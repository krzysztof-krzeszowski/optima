import copy
import gc
import numpy
import pyfits

from wx.lib import pubsub

DEFAULT_BINNING = 1.0

class Channel(object):
    def __init__(self, channel_number):
        self.channel_number = channel_number
        self.data = None
        self.second_counts = []

    def add_counts(self, second, counts):
        seconds = [v[0] for v in self.second_counts]
        if second in seconds:
            index = seconds.index(second)
            self.second_counts[index][1] += counts
        else:
            self.second_counts.append([second, counts])

    def equalise(self, factor, regions):
        if self.channel_number != '00':
            self.counts[1] *= factor / self.get_flat(regions)
        return self.counts

    def get_counts(self, binning=DEFAULT_BINNING):
        if self.counts and self.binning == binning:
            return self.counts
        self.binning = binning
        times = self.seconds + self.subseconds
        # +2 because last value doesn't get into arange and last photon has t = last second + subsecond
        bins = numpy.arange(self.seconds[0], self.seconds[-1] + 2, binning)
        y, x = numpy.histogram(times, bins)
        self.counts = [x[:-1], y]
        if self.dark_substracted:
            self.counts = self.substract_dark()
        return self.counts

    def get_flat(self, regions):
        if regions is not None:
            self.f_1 = regions[0]
            self.f_2 = regions[1]
            self.f_3 = regions[2]
            self.f_4 = regions[3]
        f_1 = int((self.f_1 - self.t_0) / self.binning)
        f_2 = int((self.f_2 + 1 - self.t_0) / self.binning)
        f_3 = int((self.f_3 - self.t_0) / self.binning)
        f_4 = int((self.f_4 + 1 - self.t_0) / self.binning)
        f1 = numpy.mean(self.counts[1][f_1:f_2])
        f2 = numpy.mean(self.counts[1][f_3:f_4])
        f = (f1 + f2) / 2
        return f

    def prepare_data(self, min_last_sec):
        self.data = numpy.array(([v[0] for v in self.second_counts if v[0] <= min_last_sec], [v[1] for v in self.second_counts if v[0] <= min_last_sec]))
        del(self.second_counts)

    def substract_dark(self, regions=None):
        if regions is not None:
            self.r_1 = regions[0]
            self.r_2 = regions[1]
            self.r_3 = regions[2]
            self.r_4 = regions[3]
        r_1 = int((self.r_1 - self.t_0) / self.binning)
        r_2 = int((self.r_2 + 1 - self.t_0) / self.binning)
        r_3 = int((self.r_3 - self.t_0) / self.binning)
        r_4 = int((self.r_4 + 1 - self.t_0) / self.binning)
        d1 = numpy.mean(self.counts[1][r_1:r_2])
        d2 = numpy.mean(self.counts[1][r_3:r_4])
        d = (d1 + d2) / 2
        self.counts[1] -= d
        self.dark_substracted = True
        return self.counts

class DataFrame(object):
    def __init__(self):
        self.files = {}
        self.channels = {}
        self.counts = []
        self.header = None

    def clear_everything(self):
        for ch in self.files.keys():
            for f in self.files[ch].values():
                f.close()
            del(self.files[ch])
        del(self.files)
        self.files = {}
        self.channels = {}

    def equalise(self, regions):
        channels = self.get_channel_numbers()
        factor = self.channels['00'].get_flat(regions)
        for channel in channels:
            self.counts[channel] = self.channels[channel].equalise(factor, regions)
        pubsub.Publisher().sendMessage('REPLOT')

    def get_binning(self):
        return self.parent.get_binning()

    def get_channel_numbers(self):
        return self.channels.keys()

    def get_counts(self, binning=DEFAULT_BINNING):
        if self.counts and binning == self.binning:
            return self.counts
        counts = {}
        binning = float(binning)
        self.binning = binning
        channels = self.get_channel_numbers()
        for channel in channels:
            counts[channel] = self.channels[channel].get_counts(self.binning)
        self.counts = counts
        return self.counts

    def get_channel_and_n(self, file_name):
        # file_name i.e.: 1rxsj211336.1_phot_ch01_2011-07-05T22-34-07_n01.fits
        # only channel number and n are needed
        a = file_name.split('/')[-1].split('.fits')[0].split('_')
        ch = '%02d' % (int(a[-3][-2:]) - 1)
        print 'ch: ', ch
        n = a[-1][-2:]
        return ch, n

    def get_dark_or_flat_value(self, channel, regions):
        """ returns dark or flat levels counts per second """
        offset = self.channels[channel].data[0][0]
        d1 = numpy.mean(self.channels[channel].data[1][regions[0] - offset:regions[1] - offset])
        d2 = numpy.mean(self.channels[channel].data[1][regions[2] - offset:regions[3] - offset])
        return (d1 + d2) / 2

    def get_data_to_save(self):
        to_return = {}
        to_return['MJD-OBS'] = self.mjd_obs
        to_return['data'] = self.get_substracted_data()
        return to_return

    def get_header_txt(self):
        header = (  ('Target',  self.header['TARGETID']),
                    ('Telescope', self.header['TELESCOP']),
                    ('Observing date [UTC]', self.header['DATE-OBS']),
                    ('MJD zero point', self.header['MJD-OBS']),
                    ('Time unit', 's'),
                    ('Barycentric correction', 'FALSE'),
                    )
        return header

    def get_lc_data(self, channels, source_channels, binning, dark_regions, flat_regions):
        """ return resulting light curve from all source channels """
        lc = []

        # get dark substracted and flat equalised data for all channels
        processed_data = self.get_processed_data(channels=channels, binning=binning, dark_regions=dark_regions, flat_regions=flat_regions)

        # prepare background channels: copy processed data and remove source channels
        background_channels = copy.deepcopy(processed_data)
        for s_c in source_channels:
            background_channels.pop(s_c)

        # get lightcurves for all the source channels
        for s_c in source_channels:
            # get source channel
            source_channel = processed_data[s_c]
            # append the difference of source channel and mean background to lc list
            lc.append(source_channel - numpy.mean(background_channels.values(), axis=0))

        # sum all lcs
        lc = numpy.sum(lc, axis=0)

        # crop the lc
        n_points = binning * (lc.__len__() // binning)
        n_points = int(n_points)
        lc = lc[:n_points]

        binning = int(binning)

        lc = [numpy.sum(lc[x:x + binning]) for x in xrange(0, n_points, binning)]

        # return the sum of all light curves
        channels_txt = '+'.join(source_channels)
        return {channels_txt: lc}

    def get_mean_for_bins(self):
        channels = self.channels
        del(channels['00'])
        channels = channels.values()
        to_return = numpy.mean([c.counts[1] for c in channels], 0)
        return to_return

    def get_range(self):
        mi = 9999999
        ma = -9999999
        for channel in self.channels.values():
            if mi > channel.data[0][0]:
                mi = channel.data[0][0]
            if ma < channel.data[1][-1] + 1:
                ma = channel.data[1][-1] + 1
        return mi, ma

    def get_processed_data(self, channels, binning, dark_regions, flat_regions, **kwargs):
        data = {}
        darks = {}
        flats = {}
        first_channel = channels[0]
        for c in channels:
            base_dark = self.get_dark_or_flat_value(c, dark_regions)
            base_flat = self.get_dark_or_flat_value(c, flat_regions) - base_dark
            darks[c] = base_dark
            flats[c] = base_flat
        first_flat = flats[first_channel]
        for c in channels:
            data[c] = (self.channels[c].data[1] - darks[c]) * first_flat / flats[c]
        return data

    def get_raw_data(self, channels, binning, **kwargs):
        data = {}
        for c in channels:
            data[c] = self.channels[c].data[1]
        return data

    def get_substracted_data(self):
        to_return = numpy.array(self.channels['00'].counts)
        to_return[1] -= self.get_mean_for_bins()
        return to_return
        
    def open_files(self, file_list):
        """ opening files clears everything (see: optima.py) """
        for f in file_list:
            ch, n = self.get_channel_and_n(f)

            print 'opening channel %s file number %s (%s)' % (ch, n, f)

            fits = pyfits.open(f)

            if not self.header:
                self.header = fits[1].header

            seconds = fits[1].data.field(0)
            first_second = int(seconds[0])
            last_second = int(seconds[-1])

            # skipping first second, very bad hack :/
            if first_second == 1:
                first_second = 2

            if not self.channels.has_key(ch):
                self.channels[ch] = Channel(ch)

            for i in xrange(first_second, last_second + 1):
                tmp = seconds[seconds == i]
                self.channels[ch].add_counts(i, tmp.__len__())

            fits.close()
            fits = pyfits.open(f)

            if not self.files.has_key(ch):
                self.files[ch] = {}
    
            self.files[ch][n] = fits

            seconds = None
            subseconds = None
            del(seconds)
            del(subseconds)
            gc.collect()

        # finished opening files
        print '\n\nReading done...\n\n'

        min_last_sec = 9e19

        # check what is the last second in each of channels and trim all channels
        for c in self.channels.values():
            if c.second_counts[-1][0] < min_last_sec:
                min_last_sec = c.second_counts[-1][0]

        min_last_sec -= 1

        for c in self.channels.values():
            c.prepare_data(min_last_sec)

    def save_data(self, name, extension, data_type, channels, source_channels, dark_regions, flat_regions, binning):
        # take data_type and invoke appropiate property to get particular data type (raw, processed, lc)
        data_to_save = getattr(self, 'get_%s_data' % data_type)(channels=channels, source_channels=source_channels, binning=binning, dark_regions=dark_regions, flat_regions=flat_regions)
        # write data_to_save to the file based on extension (txt, fits)
        getattr(self, 'write_%s' % extension)(data_to_save=data_to_save, file_name=name, binning=binning, source_channels=source_channels)

    def save_data_old(self, channels, dark_regions, flat_regions, binning):
        channels.sort()
        dark_regions.sort()
        flat_regions.sort()
        darks = {}
        flats = {}
        data = {}
        for c in channels:
            base_dark = self.get_dark_or_flat_value(c, dark_regions)
            base_flat = self.get_dark_or_flat_value(c, flat_regions) - base_dark
            darks[c] = base_dark
            flats[c] = base_flat
        first_flat = flats['00']
        first_channel = self.channels['00'].data[1] - darks['00']
        for c in channels[1:]: # without first channel
            data[c] = (self.channels[c].data[1] - darks[c]) * first_flat / flats[c]

        f = open('channels.txt', 'w')
        f.write('# ' + ' '.join(channels) + '\n')
        for i, v in enumerate(zip(*data.values())):
            f.write(str(first_channel[i]) + ' ' + ' '.join(str(i) for i in v) + '\n')
        f.close()

        print 'channels written'

        to_save = first_channel - numpy.mean(data.values(), axis=0)
        f = open('ascii.txt', 'w')
        for v in to_save:
            f.write(str(v) + '\n')
        f.close()

    def save_raw_data(self, channels, dark_regions, flat_regions, binning):
        channels.sort()

        data = {}

        for c in channels:
            data[c] = self.channels[c].data[1]

        name = self.files.values()[0].values()[0]._HDUList__file.name.split('/')[-1].split('.')[0]
        f = open(name + '.raw.txt', 'w')
        f.write('# ' + ' '.join(channels) + '\n')
        for i, v in enumerate(zip(*[data[c] for c in channels])):
            f.write(' '.join(str(i) for i in v) + '\n')
        f.close()

    def substract_dark(self, regions):
        channels = self.get_channel_numbers()
        for channel in channels:
            self.counts[channel] = self.channels[channel].substract_dark(regions)
        pubsub.Publisher().sendMessage('REPLOT')

    def write_fits(self, data_to_save, file_name, binning, source_channels):
        header_keys_to_rewrite = ('TARGETID', 'OBSERVER', 'COMMENT', 'TYPE', 'OBSMODE', 'TELESCOP',
                                    'DATE', 'RA', 'DEC', 'EQUINOX', 'INSTRUME', 'EXPTIME', 'MJD-OBS',
                                    'DATE-OBS', 'UTCSOURC', 'AIRM', 'REG-0X01', 'REG-0X02', 'REG-0X03',
                                    'REG-0X04', 'REG-0X05', 'REG-0X06', 'REG-0X07', 'REG-0X08', 'REG-0X09',
                                    'REG-0X0A', 'REG-0X0B', 'REG-0X0C', 'REG-0X0D', 'REG-0X0E', 'LATITUDE',
                                    'LONGITUD', 'ALTITUDE', 'CALPOLAN',)
        channels = data_to_save.keys()
        channels.sort()
        columns = []
        seconds = self.channels.values()[0].data[0]
        columns.append(pyfits.Column(name='TIME', format='E', array=seconds))
        for i, ch in enumerate(channels):
            if i == 0:
                c = pyfits.Column(name='RATE', format='E', array=data_to_save[ch])
            else:
                c = pyfits.Column(name='ch' + ch, format='E', array=data_to_save[ch])
            columns.append(c)

        coldefs = pyfits.ColDefs([c for c in columns])
        tbhdu = pyfits.new_table(coldefs)

        for key in header_keys_to_rewrite:
            tbhdu.header.update(key, self.header[key])

        # keywords required by xronos
        # http://heasarc.gsfc.nasa.gov/docs/xanadu/xronos/manual/node5.html

        tbhdu.header.update('MJDREF', tbhdu.header['MJD-OBS'])
        tbhdu.header.update('TSTART', seconds[0])
        tbhdu.header.update('TSTOP', seconds[-1])
        tbhdu.header.update('TIMEZERO', 0)
        tbhdu.header.update('TIMESYS', tbhdu.header['MJD-OBS'])
        tbhdu.header.update('TIMEUNIT', 's')
        tbhdu.header.update('TIMEDEL', seconds[1] - seconds[0])
        
        channels_txt = '+'.join(source_channels)
        tbhdu.header.update('CHANNELS', channels_txt)

        hdu = pyfits.PrimaryHDU()
        thdulist = pyfits.HDUList([hdu, tbhdu])

        tbhdu.writeto(file_name, clobber=True)

    def write_txt(self, data_to_save, file_name, binning, **kwargs):
        seconds = self.channels.values()[0].data[0]
        seconds = range(seconds[0], seconds[0] + seconds.__len__(), int(binning))
        header = self.get_header_txt()
        channels = data_to_save.keys()
        channels.sort()
        f = open(file_name, 'w')
        for h in header:
            f.write('# ' + h[0] + ': ' + str(h[1]) + '\n')
        f.write('# Binning: ' + str(binning) + '\n')
        f.write('# Channels: ' + ' '.join(channels) + '\n')
        for s, vals in zip(seconds, zip(*[data_to_save[c] for c in channels])):
            f.write(str(s + 0.5 * binning) + ' ' + ' '.join([str(v) for v in vals]) + '\n')
        f.close()
