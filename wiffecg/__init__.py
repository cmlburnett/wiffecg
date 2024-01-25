import datetime

import wiff
import matplotlib
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages

import pyiworxekgedfimport

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
					print([f[0], page,r, len(times), len(vals), len(vals[0]), len(vals[1])])
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

