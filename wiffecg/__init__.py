import wiff
import matplotlib
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages

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

	def ExportToPDF(self, fname, speed=200):
		"""
		Export EKG tracing to a PDF file at @fname.
		Save at @speed mm/sec (default 25).
		"""

		leads = self.GetLeads()
		width = 8.0

		freq = self.wiff.recording[1].sampling
		samps = (width * 25.4) / speed * freq
		step = int(samps)

		pyplot.rcParams["figure.figsize"] = [width, 10.0]
		pyplot.rcParams["figure.autolayout"] = True

		with PdfPages(fname) as pp:
			fig,axes = pyplot.subplots(len(leads))

			times = []
			vals = [[] for _ in range(len(leads))]
			for f in self.wiff.recording[1].GetAllFrames():
				page,r = divmod(f[0]+1, step)
				if r == 0:
					print([f[0], page,r, len(times), len(vals), len(vals[0]), len(vals[1])])
					for i in range(len(leads)):
						lead = leads[i]
						axes[i].cla()
						axes[i].set_ylabel(lead)
						axes[i].plot(times, vals[i])
						axes[i].xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
						axes[i].xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("{x:.2f}"))
						axes[i].xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.025))
					pp.savefig(fig)
					pyplot.cla()
					for v in vals:
						v.clear()
					times.clear()

				times.append(f[0]/freq)
				for i in range(len(leads)):
					vals[i].append(f[2][i])

