import datetime
import enum
import os
import pickle
import zipfile

import wiff
import matplotlib
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages

import pyiworxekgedfimport
from pyzestyecg import pyzestyecg

class WIFFECG:
	def __init__(self, fname):
		self._fname = fname
		self._wiff = wiff.open(fname)

	def close(self):
		self._wiff.close()

	@property
	def Filename(self):
		return self._fname

	@property
	def wiff(self):
		return self._wiff

	def __enter__(self):
		return self

	def __exit__(self, *args):
		self.close()
		return False

	def GetDuration(self):
		delta = self.wiff.recording[1].end - self.wiff.recording[1].start
		return delta.total_seconds()

	def Validate(self):
		"""
		WIFF must meet certain conditions to be considered an ECG.
		"""

		leads = self.GetLeads()
		for lead in leads:
			if lead.startswith("Lead "):
				lead = lead.split('Lead ')[1]
			if lead not in ('I', 'II', 'III', 'aVR', 'aVL', 'aVF'):
				raise ValueError("Lead '%s' not acceptable for ECG" % lead)

		if len(self.wiff.recording) != 1:
			raise ValueError("Multiple recordings (%d) found in a single file, not supported" % len(self.wiff.recording))

		ms = self.wiff.meta.find(None, 'Recording.Type')
		if not len(ms):
			raise ValueError("No meta value for Recording.Type")
		if not ms[0].value.startswith('EKG'):
			raise ValueError("Recroding.Type meta ('%s') does not start with 'EKG'" % ms[0].value)

	def GetLeads(self):
		"""
		Get all the leads (channels) present.
		"""

		ret = []
		for chan in self.wiff.channel.values():
			ret.append(chan.name)

		return ret

	@staticmethod
	def ImportFromIWorxEDF(inpath, outpath):
		"""
		Using IWorx's LabScribe, export an EKG recording as an EDF file.
		This uses pyiworxekgedfimport to import the data into a WIFF file.
		@inpath is file path to the EDF file.
		@outpath is the file path to save the WIFF file at.
		"""

		with pyiworxekgedfimport.EDFReader.open(inpath) as o:
			props = {
				'start': o.Start,
				'end': o.Start + datetime.timedelta(seconds=o.Duration),
				'description': 'iWorx LabScribe EDF file',
				'fs': 2000, # FIXME
				'channels': [],
			}

			for idx,ch in enumerate(o.Signals):
				if ch['Label'] == 'EDF Annotations':
					continue

				props['channels'].append({
					'idx': idx,
					'name': ch['Label'],
					'bits': 16,
					'unit': ch['PhysicalDimension'],
					'digitalminvalue': ch['DigitalMinimum'],
					'digitalmaxvalue': ch['DigitalMaximum'],
					'analogminvalue': ch['PhysicalMinimum'],
					'analogmaxvalue': ch['PhysicalMaximum'],
					'comment': 'Physical (%d,%d) to (%d,%d)' % (ch['PhysicalMinimum'],ch['PhysicalMaximum'], ch['DigitalMinimum'],ch['DigitalMaximum']),
				})

			# Save data to file
			o.writeWIFF(outpath, props)

	def ExportToPDF(self, fname, speed=100):
		"""
		Export EKG tracing to a PDF file at @fname.
		Save at @speed mm/sec (default 25).
		"""

		leads = self.GetLeads()
		width = 8.0

		# Determine how many samples per page
		freq = self.wiff.recording[1].sampling
		samps = (width * 25.4) / speed * freq
		step = int(samps)

		pyplot.rcParams["figure.figsize"] = [width, 10.0]
		pyplot.rcParams["figure.autolayout"] = True
		pyplot.rcParams["xtick.labelsize"] = 'small'

		with PdfPages(fname) as pp:
			fig,axes = pyplot.subplots(len(leads))

			times = []
			vals = [[] for _ in range(len(leads))]
			for f in self.wiff.recording[1].GetAllFrames():
				# Every time the remainder trips zero, that's a new page so save it

				# First frame is zero, but no need to save it yet so skirt around this by adding one
				page,r = divmod(f[0]+1, step)

				# Save new page of data
				if r == 0:
					for i in range(len(leads)):
						lead = leads[i]
						if lead.startswith("Lead "):
							lead = lead.split("Lead ")[1]
						axes[i].cla()
						axes[i].set_ylabel(lead)
						axes[i].plot(times, vals[i])

						# TODO: These should change depending on the speed
						axes[i].xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.5))
						axes[i].xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter("{x:.1f}"))
						axes[i].xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.05))

						# Only show label on bottom subplot
						if i+1 == len(leads):
							axes[i].set_xlabel("Time (sec)")

					# Save the figure, clear it, and clear the data (rather than making a new figure each time == memory wasteful)
					pp.savefig(fig)
					pyplot.cla()
					for v in vals:
						v.clear()
					times.clear()

				# Append time and data
				times.append(f[0]/freq)
				for i in range(len(leads)):
					vals[i].append(f[2][i])

	def ProcessRRToZip(self, fname, params, intervals, savepng=True, savepdf=True):
		"""
		Process the RR intervals and use a zip file to store all the results.
		@intervals should be a dictionary
			'user'		user defined intervals, a dictionary mapping to a list of (start,stop) tuples
			'noise'		list of (start,stop) tuples that indicate noise
			'ignore'	list of (start,stop) tuples that are ignored entirely
		@savepng indicates if PNG images of the EKG should be saved
		@savepdf indicates if PDF file of the EKG with peaks should be saved
		"""

		# Get sampling rate for the recording
		samp = self.wiff.recording[1].sampling

		# intervals are in seconds, convert to framex by multiplying by the sampling rate
		# Round to the nearest frame
		intervals_user_frames = {}
		for k,v in intervals['user'].items():
			intervals_user_frames[k] = [(round(_[0]*samp),round(_[1]*samp)) for _ in v]

		intervals_ignore_frames = []
		for v in intervals['ignore']:
			intervals_ignore_frames.append( (round(v[0]*samp),round(v[1]*samp)) )

		intervals_noise_frames = []
		for v in intervals['noise']:
			intervals_noise_frames.append( (round(v[0]*samp),round(v[1]*samp)) )

		with ZipMan(fname) as z:
			if z.IsStateEmpty:
				# Start processing from the start, do any initializing
				z.State['Channels'] = self.GetLeads()
				z.State['Potentials'] = None
				z.State['Peaks'] = None
				z.State['Correlate'] = None
				z.State['Keep'] = None
				z.State['Remove'] = None
				z.State['UserFilter'] = None
				z.SetStatInitialized()
				z.SaveState()

			elif z.IsStateInitialized:
				p = pyzestyecg(self.wiff, params)
				potentials, peaks = p.GetPotentialsAndPeaks(intervals_ignore_frames, intervals_noise_frames)

				z.State['Potentials'] = potentials
				z.State['Peaks'] = peaks
				z.SetStateCorrelate()
				z.SaveState()

			elif z.IsStateCorrelate:
				chans = z.State['Channels']
				potentials = z.State['Potentials']
				peaks = z.State['Peaks']

				p = pyzestyecg(self.wiff)
				correlate = p.GetCorrelate(chan, potentials, peaks)

				z.State['Correlate'] = correlate
				z.SetStateKeepKeys()
				z.SaveState()

			elif z.IsStateKeepKeys:
				chans = z.State['Channels']
				peaks = z.State['Peaks']

				p = pyzestyecg(self.wiff)
				points, keep = p.GetKeepKeys(chans, peaks)

				z.State['Keep'] = keep
				z.State['Points'] = points
				z.SetStateRemoveKeys()
				z.SaveState()

			elif z.IsStateRemoveKeys:
				chans = z.State['Channels']
				correlate = z.State['Correlate']
				points = z.State['Points']
				keep = z.State['Keep']

				p = pyzestyecg(self.wiff)
				remove = p.GetRemoveKeys(chans, correlate, points, keep)

				z.State['Remove'] = remove
				z.SetStateUserFilter()
				z.SaveState()

			elif z.IsStateUserFilter:
				chans = z.State['Channels']
				points = z.State['Points']
				keep = z.State['Keep']
				remove = z.State['Remove']
				user = z.State['UserFilter']

				p = pyzestyecg(self.wiff)
				p.CheckUserFilter(chans, points, keep, remove, user)

				z.SetStateCalculateRR()
				z.SaveState()

			elif z.IsStateCalculateRR:
				chans = z.State['Channels']
				correlate = z.State['Correlate']
				keep = z.State['Keep']
				remove = z.State['Remove']
				user = z.State['UserFilter']

				p = pyzestyecg(self.wiff)
				p.CalculateRR(chans, keep, remove, user, intervals_user_frames, intervals_noise_frames, params)

				if savepng:
					z.SetStateSavePNG()
				elif savepdf:
					z.SetStateSavePDF()
				else:
					z.SetStateCompleted()
				z.SaveState()


			elif z.IsStateSavePNG:
				print("Save as PNG")

				if savepdf:
					z.SetStateSavePDF()
				else:
					z.SetStateCompleted()
				z.SaveState()

			elif z.IsStateSavePDF:
				print("Save as PDF")

				z.SetStateCompleted()
				z.SaveState()

			elif z.IsStateCompleted:
				raise NotImplementedError("Should not reach this point")

			else:
				raise NotImplementedError("Should not reach this point")

		return z.IsStateCompleted

class ZipMan:
	"""
	ZipMan -- Zip Manager
	Wraps a zip file to store state in a pickle file and other output files.
	"""

	class ProcessingStateEnum(enum.IntEnum):
		EMPTY = 0
		INITIALIZED = 1
		CORRELATE = 2
		KEEPKEYS = 3
		REMOVEKEYS = 4
		USERFILTER = 5
		CALCULATERR = 6
		SAVEPNG = 7
		SAVEPDF = 8
		COMPLETED = 100

		ERROR = 1000

	def __init__(self, fname):
		self._fname = fname
		self._zip = None
		self._state = {
			'state': ZipMan.ProcessingStateEnum.EMPTY,
			'error': None,
		}

	@property
	def Filename(self):
		return self._fname

	@property
	def Zip(self):
		return self._zip

	@property
	def State(self):
		return self._state


	@property
	def IsStateEmpty(self):
		return self._state['state'] == ZipMan.ProcessingStateEnum.EMPTY

	@property
	def IsStateInitialized(self):
		return self._state['state'] == ZipMan.ProcessingStateEnum.INITIALIZED

	@property
	def IsStateCorrelate(self):
		return self._state['state'] == ZipMan.ProcessingStateEnum.CORRELATE

	@property
	def IsStateKeepKeys(self):
		return self._state['state'] == ZipMan.ProcessingStateEnum.KEEPKEYS

	@property
	def IsStateRemoveKeys(self):
		return self._state['state'] == ZipMan.ProcessingStateEnum.REMOVEKEYS

	@property
	def IsStateUserFilter(self):
		return self._state['state'] == ZipMan.ProcessingStateEnum.USERFILTER

	@property
	def IsStateCalculateRR(self):
		return self._state['state'] == ZipMan.ProcessingStateEnum.CALCULATERR

	@property
	def IsStateSavePNG(self):
		return self._state['state'] == ZipMan.ProcessingStateEnum.SAVEPNG

	@property
	def IsStateSavePDF(self):
		return self._state['state'] == ZipMan.ProcessingStateEnum.SAVEPDF

	@property
	def IsStateCompleted(self):
		return self._state['state'] == ZipMan.ProcessingStateEnum.COMPLETED

	@property
	def IsStateNonError(self):
		return self._state['state'] != ZipMan.ProcessingStateEnum.ERROR

	@property
	def IsStateError(self):
		return self._state['state'] == ZipMan.ProcessingStateEnum.ERROR


	def SetStatEmpty(self):
		self._state['state'] = ZipMan.ProcessingStateEnum.EMPTY

	def SetStatInitialized(self):
		self._state['state'] = ZipMan.ProcessingStateEnum.INITIALIZED

	def SetStateCorrelate(self):
		self._state['state'] = ZipMan.ProcessingStateEnum.CORRELATE

	def SetStateKeepKeys(self):
		self._state['state'] = ZipMan.ProcessingStateEnum.KEEPKEYS

	def SetStateRemoveKeys(self):
		self._state['state'] = ZipMan.ProcessingStateEnum.REMOVEKEYS

	def SetStateUserFilter(self):
		self._state['state'] = ZipMan.ProcessingStateEnum.USERFILTER

	def SetStateCalculateRR(self):
		self._state['state'] = ZipMan.ProcessingStateEnum.CALCULATERR

	def SetStateSavePNG(self):
		self._state['state'] = ZipMan.ProcessingStateEnum.SAVEPNG

	def SetStateSavePDF(self):
		self._state['state'] = ZipMan.ProcessingStateEnum.SAVEPDF

	def SetStateCompleted(self):
		self._state['state'] = ZipMan.ProcessingStateEnum.COMPLETED

	def SetStateError(self):
		self._state['state'] = ZipMan.ProcessingStateEnum.ERROR

	def SetError(self, e, tb, msg, dat=None):
		"""
		Set the state to error and record the exception instance @e, the traceback @tb, and user message @msg.
		Optional is a data object @dat to include (eg, a dictionary of variables)
		"""
		self.SetStateError()
		self._state['error'] = {'exception': e, 'traceback': tb, 'message': msg, 'data': dat}

	def open(self):
		if self._zip:
			raise Exception("File is already open")

		self._zip = zipfile.ZipFile(self.Filename, 'a')

		# Write empty state if not existant or read what's there
		if 'state.pypickle' not in self.Zip.namelist():
			with self.Zip.open('state.pypickle', 'w') as f:
				pickle.dump(self.State, f)
		else:
			with self.Zip.open('state.pypickle', 'r') as f:
				self._state = pickle.load(f)

	def close(self):
		if self._zip is None:
			raise Exception("File is not open, cannot close it")

		self._zip.close()
		self._zip = None

	def SaveState(self):
		# Save the state
		with self.Zip.open('state.pypickle', 'w') as f:
			pickle.dump(self.State, f)

	def __enter__(self):
		self.open()

	def __exit__(self, *kargs):
		if self.Zip:
			self.close()
		return False

