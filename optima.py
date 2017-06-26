#!/usr/bin/env python

import matplotlib
matplotlib.use('WXAgg')
import numpy
import os
import pylab
import wx

from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar

from new_data import DataFrame

COLORS = (
            wx.Colour(254, 0, 0), # ch00
            wx.Colour(252, 127, 128), # ch01
            wx.Colour(0, 253, 0), # ch02
            wx.Colour(128, 253, 128), # ch03
            wx.Colour(253, 253, 0), # ch04
            wx.Colour(253, 253, 127), # ch05
            wx.Colour(127, 126, 252), # ch06
            wx.Colour(180, 180, 253), # ch07
            wx.Colour(199, 199, 199), # ch08
            wx.Colour(253, 169, 0), # ch09
            wx.Colour(84, 168, 252), # ch10
            wx.Colour(252, 126, 253),# ch11
        )
DEFAULT_BINNING = 1
DEFAULT_PLOT_RANGE = 600
PLOT_OFFSET = 50
Y_MAX_SCALLING = 1.2

DARK_POSITION_1_OFFSET = 10
DARK_POSITION_2_OFFSET = 100
DARK_POSITION_3_OFFSET = -100
DARK_POSITION_4_OFFSET = -10
FLAT_POSITION_1_OFFSET = 110
FLAT_POSITION_2_OFFSET = 200
FLAT_POSITION_3_OFFSET = -200
FLAT_POSITION_4_OFFSET = -110

# offsets for moving dark and flat lines and screen scroll
DARK_LINE_OFFSET = 100
FLAT_LINE_OFFSET = 100

class ControlPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)

        self.parent = parent

        self.colour_pickers = {}
        self.checkboxes = {}
        self.source_checkboxes = {}
        self.darks = []
        self.flats = []

        self.working_path = None
        self.working_dir = None

        sizer = wx.BoxSizer(wx.HORIZONTAL)

        # left sizer
        left_panel = wx.BoxSizer(wx.VERTICAL)

        open_files_button = wx.Button(self, wx.ID_ANY, 'Open files')
        open_files_button.Bind(wx.EVT_BUTTON, self.on_open_files)
        left_panel.Add(open_files_button, 0, wx.EXPAND)

        left_panel.Add(wx.StaticText(self, wx.ID_ANY, 'Binning [s]'), 0, wx.EXPAND)
        self.binning = wx.TextCtrl(self, wx.ID_ANY, '%i' % DEFAULT_BINNING, style=wx.TE_RIGHT)
        self.binning.Enable(0)
        left_panel.Add(self.binning, 0, wx.EXPAND)

        left_panel.Add(wx.StaticText(self, wx.ID_ANY, 'Save data'), 0, wx.ALIGN_CENTER)

        save_sizer = wx.BoxSizer(wx.HORIZONTAL)

        save_ascii_sizer = wx.BoxSizer(wx.VERTICAL)

        save_ascii_sizer.Add(wx.StaticText(self, wx.ID_ANY, 'ASCII'), 0, wx.ALIGN_CENTER)
        self.save_raw_ascii_button = wx.Button(self, wx.ID_ANY, 'Raw', name='txt_raw')
        self.save_raw_ascii_button.Bind(wx.EVT_BUTTON, self.on_save)
        self.save_raw_ascii_button.Enable(0)
        save_ascii_sizer.Add(self.save_raw_ascii_button, 0, wx.EXPAND)

        self.save_processed_ascii_button = wx.Button(self, wx.ID_ANY, 'Processed', name='txt_processed')
        self.save_processed_ascii_button.Bind(wx.EVT_BUTTON, self.on_save)
        self.save_processed_ascii_button.Enable(0)
        save_ascii_sizer.Add(self.save_processed_ascii_button, 0, wx.EXPAND)

        self.save_lc_ascii_button = wx.Button(self, wx.ID_ANY, 'Light curve', name='txt_lc')
        self.save_lc_ascii_button.Bind(wx.EVT_BUTTON, self.on_save)
        self.save_lc_ascii_button.Enable(0)
        save_ascii_sizer.Add(self.save_lc_ascii_button, 0, wx.EXPAND)

        save_sizer.Add(save_ascii_sizer, 0, wx.EXPAND)

        save_fits_sizer = wx.BoxSizer(wx.VERTICAL)

        save_fits_sizer.Add(wx.StaticText(self, wx.ID_ANY, 'FITS'), 0, wx.ALIGN_CENTER)
        self.save_raw_fits_button = wx.Button(self, wx.ID_ANY, 'Raw', name='fits_raw')
        self.save_raw_fits_button.Bind(wx.EVT_BUTTON, self.on_save)
        self.save_raw_fits_button.Enable(0)
        save_fits_sizer.Add(self.save_raw_fits_button, 0, wx.EXPAND)

        self.save_processed_fits_button = wx.Button(self, wx.ID_ANY, 'Processed', name='fits_processed')
        self.save_processed_fits_button.Bind(wx.EVT_BUTTON, self.on_save)
        self.save_processed_fits_button.Enable(0)
        save_fits_sizer.Add(self.save_processed_fits_button, 0, wx.EXPAND)

        self.save_lc_fits_button = wx.Button(self, wx.ID_ANY, 'Light curve', name='fits_lc')
        self.save_lc_fits_button.Bind(wx.EVT_BUTTON, self.on_save)
        self.save_lc_fits_button.Enable(0)
        save_fits_sizer.Add(self.save_lc_fits_button, 0, wx.EXPAND)

        save_sizer.Add(save_fits_sizer, 0, wx.EXPAND)
        left_panel.Add(save_sizer, 0, wx.EXPAND)

        self.plot_lc_button = wx.Button(self, wx.ID_ANY, 'Plot light curve')
        self.plot_lc_button.Bind(wx.EVT_BUTTON, self.on_plot_lc)
        self.plot_lc_button.Enable(0)
        left_panel.Add(self.plot_lc_button, 0, wx.EXPAND)

        self.subsecond_binning = wx.Button(self, wx.ID_ANY, 'Subsecond binning')
        self.subsecond_binning.Bind(wx.EVT_BUTTON, self.on_subsecond_binning)
        self.subsecond_binning.Enable(0)
        left_panel.Add(self.subsecond_binning, 0, wx.EXPAND)

        # channel panel
        channel_panel = wx.GridSizer(rows=4, cols=3, vgap=5, hgap=5)

        for i in xrange(12):
            str_ch_no = '%02d' % i
            channel_sizer = wx.BoxSizer(wx.HORIZONTAL)
            static_text = wx.StaticText(self, wx.ID_ANY, 'Ch #' + str_ch_no)
            colour_picker = wx.ColourPickerCtrl(self, wx.ID_ANY, COLORS[i], name='colour_ch_%s' % str_ch_no)
            self.Bind(wx.EVT_COLOURPICKER_CHANGED, self.on_colour_picker_change, colour_picker)
            self.colour_pickers[str_ch_no] = colour_picker
            colour_picker.Enable(0)
            checkbox = wx.CheckBox(self, wx.ID_ANY, name='plot_ch_%s' % str_ch_no, label='Plot')
            self.Bind(wx.EVT_CHECKBOX, self.on_checkbox_change, checkbox)
            self.checkboxes[str_ch_no] = checkbox
            checkbox.Enable(0)
            channel_sizer.Add(static_text, 0)
            channel_sizer.Add(colour_picker, 0)
            channel_sizer.Add(checkbox, 0)
            channel_panel.Add(channel_sizer, 0)

        # source panel
        static_box = wx.StaticBox(self, -1, 'Source in:')
        static_box_sizer = wx.StaticBoxSizer(static_box, wx.VERTICAL)
        source_panel = wx.GridSizer(rows=4, cols=3, vgap=20)
        for i in xrange(12):
            str_ch_no = '%02d' % i
            source_channel_sizer = wx.BoxSizer(wx.HORIZONTAL)
            checkbox = wx.CheckBox(self, wx.ID_ANY, name='source_in_ch_%s' % str_ch_no, label='Ch #' + str_ch_no)
            self.Bind(wx.EVT_CHECKBOX, self.on_source_checkbox_change, checkbox)
            self.source_checkboxes[str_ch_no] = checkbox
            checkbox.Enable(0)
            source_channel_sizer.Add(checkbox, 0)
            source_panel.Add(source_channel_sizer, 0)
        static_box_sizer.Add(source_panel)

        # regions panel
        regions_panel = wx.BoxSizer(wx.HORIZONTAL)

        dark_panel = wx.BoxSizer(wx.VERTICAL)
        dark_panel.Add(wx.StaticText(self, wx.ID_ANY, 'Dark regions'), 0)

        for i in xrange(4):
            dark = wx.SpinCtrl(self, wx.ID_ANY, name='dark_%s' % str(i + 1), initial=0, min=0, max=0)  
            dark.Enable(0)
            self.Bind(wx.EVT_SPINCTRL, self.on_change_dark, dark)
            self.darks.append(dark)
            dark_panel.Add(dark, 0)

        regions_panel.Add(dark_panel, 0)
    
        flat_panel = wx.BoxSizer(wx.VERTICAL)
        flat_panel.Add(wx.StaticText(self, wx.ID_ANY, 'Flat regions'), 0)

        for i in xrange(4):
            flat = wx.SpinCtrl(self, wx.ID_ANY, name='flat_%s' % str(i + 1), initial=0, min=0, max=0)  
            flat.Enable(0)
            self.Bind(wx.EVT_SPINCTRL, self.on_change_flat, flat)
            self.flats.append(flat)
            flat_panel.Add(flat, 0)

        regions_panel.AddSpacer(15)
        regions_panel.Add(flat_panel, 0)

        # add all sizers

        sizer.Add(left_panel, 0, wx.EXPAND)   
        sizer.AddSpacer(15)
        sizer.Add(channel_panel, 0, wx.EXPAND)   
        sizer.AddSpacer(15)
        sizer.Add(static_box_sizer, 0, wx.EXPAND)   
        sizer.AddSpacer(15)
        sizer.Add(regions_panel, 0, wx.EXPAND)   
        self.SetSizer(sizer)

    def clear_everything(self):
        for c in self.checkboxes.values():
            c.SetValue(0)
            c.Enable(0)
        for c in self.source_checkboxes.values():
            c.SetValue(0)
            c.Enable(0)
        for c in self.colour_pickers.values():
            c.Enable(0)
        for d in self.darks:
            d.SetRange(0, 0)
            d.SetValue(0)
            d.Enable(0)
        for f in self.flats:
            f.SetRange(0, 0)
            f.SetValue(0)
            f.Enable(0)
        self.working_path = None
        self.working_dir = None
        self.save_raw_ascii_button.Enable(0)
        self.save_processed_ascii_button.Enable(0)
        self.save_lc_ascii_button.Enable(0)
        self.save_raw_fits_button.Enable(0)
        self.save_processed_fits_button.Enable(0)
        self.save_lc_fits_button.Enable(0)
        self.plot_lc_button.Enable(0)
        self.subsecond_binning.Enable(0)

    def enable(self, channels):
        for channel in channels:
            self.checkboxes[channel].Enable(1)
            self.checkboxes[channel].SetValue(1)
            self.colour_pickers[channel].Enable(1)
            self.source_checkboxes[channel].Enable(1)
        data = self.parent.data_frame.channels.values()[0].data[0]
        start, end = data[0], data[-1]
        for dark in self.darks:
            dark.SetRange(start, end)
            dark.Enable(1)
        self.darks[0].SetValue(start + DARK_POSITION_1_OFFSET)
        self.darks[1].SetValue(start + DARK_POSITION_2_OFFSET)
        self.darks[2].SetValue(end + DARK_POSITION_3_OFFSET)
        self.darks[3].SetValue(end + DARK_POSITION_4_OFFSET)
        for flat in self.flats:
            flat.SetRange(start, end)
            flat.SetValue(start)
            flat.Enable(1)
        self.flats[0].SetValue(start + FLAT_POSITION_1_OFFSET)
        self.flats[1].SetValue(start + FLAT_POSITION_2_OFFSET)
        self.flats[2].SetValue(end + FLAT_POSITION_3_OFFSET)
        self.flats[3].SetValue(end + FLAT_POSITION_4_OFFSET)
        self.save_raw_ascii_button.Enable(1)
        self.save_processed_ascii_button.Enable(1)
        self.save_lc_ascii_button.Enable(1)
        self.save_raw_fits_button.Enable(1)
        self.save_processed_fits_button.Enable(1)
        self.save_lc_fits_button.Enable(1)
        self.plot_lc_button.Enable(1)
        self.binning.Enable(1)
        self.subsecond_binning.Enable(1)

    def get_binning(self):
        return float(self.binning.GetValue())

    def get_current_channels(self):
        channels = []
        for c in self.checkboxes.values():
            if c.GetValue():
                channels.append(str(c.GetName()).split('_')[-1])
        channels.sort()
        return channels

    def get_dark_values(self):
        darks = []
        for d in self.darks:
            darks.append(d.GetValue())
        darks.sort()
        return darks

    def get_flat_values(self):
        flats = []
        for f in self.flats:
            flats.append(f.GetValue())
        flats.sort()
        return flats

    def get_source_channels(self):
        source_channels = []
        for c in self.source_checkboxes.values():
            if c.GetValue():
                source_channels.append(str(c.GetName()).split('_')[-1])
        source_channels.sort()
        return source_channels

    def on_change_dark(self, event):
        spin = event.GetEventObject()
        n = int(str(spin.GetName()).split('_')[-1]) - 1
        v = spin.GetValue()
        self.parent.set_dark_line(n, v)

    def on_change_flat(self, event):
        spin = event.GetEventObject()
        n = int(str(spin.GetName()).split('_')[-1]) - 1
        v = spin.GetValue()
        self.parent.set_flat_line(n, v)

    def on_checkbox_change(self, event):
        checkbox = event.GetEventObject()
        name = checkbox.GetName()
        channel = str(name).split('_')[-1]
        visible = checkbox.GetValue()
        self.parent.toggle_line(channel, visible)
        if not visible:
            if self.source_checkboxes[channel].GetValue():
                self.source_checkboxes[channel].SetValue(False)

    def on_colour_picker_change(self, event):
        colour_picker = event.GetEventObject()
        name = colour_picker.GetName()
        channel = str(name).split('_')[-1]
        colours = colour_picker.GetColour()
        self.parent.toggle_line_colour(channel, colours)
        self.parent.mini_map.toggle_line_colour(channel, colours)
        
    def on_open_files(self, event):
        if not self.working_path:
            self.working_path = os.getcwd()
        dialog = OpenFilesDialog(default_dir=self.working_path)
        dialog.Center()
        if dialog.ShowModal() == wx.ID_OK:
            self.clear_everything()
            self.file_list = dialog.GetPaths()
            self.parent.open_files(file_list=dialog.GetPaths())
            self.working_path = dialog.GetDirectory()
            self.working_dir = self.working_path.split('/')[-1]
            self.parent.SetTitle('OPTIMA: ' + self.working_dir)
        dialog.Destroy()

    def on_plot_lc(self, event):
        if not self.sources_selected():
            return
        channels=self.get_current_channels()
        source_channels=self.get_source_channels()
        binning=self.get_binning()
        dark_regions = self.get_dark_values()
        flat_regions = self.get_flat_values()
        lc = self.parent.data_frame.get_lc_data(channels=channels, source_channels=source_channels, binning=binning, dark_regions=dark_regions, flat_regions=flat_regions).values()[0]
        dialog = PlotLC(lc=lc)
        dialog.Center()
        dialog.ShowModal()
        dialog.Destroy()

    def on_save(self, event):
        if not self.sources_selected():
            return
        button = event.GetEventObject()
        name = button.GetName()
        extension, data_type = name.split('_')
        if not self.working_path:
            self.working_path = os.getcwd()
        dialog = SaveFileDialog(default_dir=self.working_path, default_file=self.working_dir, data_type=data_type, extension=extension)
        dialog.Center()
        if dialog.ShowModal() == wx.ID_OK:
            channels = self.get_current_channels()
            darks = self.get_dark_values()
            flats = self.get_flat_values()
            binning = self.get_binning()
            name = dialog.GetPath()
            self.working_path = os.path.dirname(name)
            source_channels = self.get_source_channels()
            self.parent.data_frame.save_data(name=name, extension=extension, data_type=data_type, channels=channels, source_channels=source_channels, dark_regions=darks, flat_regions=flats, binning=binning)
        dialog.Destroy()

    def on_save_data(self, event):
        channels = []
        for c in self.checkboxes.values():
            if c.GetValue():
                channels.append(str(c.GetName()).split('_')[-1])
        darks = []
        for d in self.darks:
            darks.append(d.GetValue())
        flats = []
        for f in self.flats:
            flats.append(f.GetValue())
        binning = float(self.binning.GetValue())
        self.parent.data_frame.save_data(channels=channels, dark_regions=darks, flat_regions=flats, binning=binning)

    def on_save_fits(self, event):
        if not self.working_path:
            self.working_path = os.getcwd()
        dialog = SaveFileDialog(default_dir=self.working_path, default_file=self.working_dir, extension='fits')
        dialog.Center()
        if dialog.ShowModal() == wx.ID_OK:
            channels = []
            for c in self.checkboxes.values():
                if c.GetValue():
                    channels.append(str(c.GetName()).split('_')[-1])
            darks = []
            for d in self.darks:
                darks.append(d.GetValue())
            flats = []
            for f in self.flats:
                flats.append(f.GetValue())
            binning = float(self.binning.GetValue())
            name = dialog.GetPaths()
            self.parent.data_frame.save_fits(name=name, channels=channels, dark_regions=darks, flat_regions=flats, binning=binning)
        dialog.Destroy()

    def on_save_raw_data(self, event):
        channels = []
        for c in self.checkboxes.values():
            if c.GetValue():
                channels.append(str(c.GetName()).split('_')[-1])
        darks = []
        for d in self.darks:
            darks.append(d.GetValue())
        flats = []
        for f in self.flats:
            flats.append(f.GetValue())
        binning = float(self.binning.GetValue())
        self.parent.data_frame.save_raw_data(channels=channels, dark_regions=darks, flat_regions=flats, binning=binning)

    def on_source_checkbox_change(self, event):
        checkbox = event.GetEventObject()
        if checkbox.GetValue():
            name = checkbox.GetName()
            channel = str(name).split('_')[-1]
            if not self.checkboxes[channel].GetValue():
                self.checkboxes[channel].SetValue(True)
                self.parent.toggle_line(channel)

    def on_subsecond_binning(self, evnet):
        if not self.sources_selected():
            return
        if not self.working_path:
            self.working_path = os.getcwd()
        dialog = SaveFileDialog(default_dir=self.working_path, default_file=self.working_dir, data_type="sub", extension="txt")
        dialog.Center()
        if dialog.ShowModal() == wx.ID_OK:
            name = dialog.GetPath()
            f = open(name, "w")
            for file_name in self.file_list:
                f.write("file: " + file_name + "\n")
            f.close()
        dialog.Destroy()

    def sources_selected(self):
        if not self.get_source_channels():
            alert = wx.MessageDialog(self, 'No source channel selected', 'Error', style=wx.OK|wx.CENTRE)
            alert.Center()
            alert.ShowModal()
            alert.Destroy()
            return False
        else:
            return True

class MiniMapPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)

        self.parent = parent

        self.dpi = 100
        self.fig = Figure((1.0, 1.0), dpi=self.dpi)

        self.fig.subplots_adjust(left=0.05, right=0.97, bottom=0.16, top=0.93)

        self.axes = self.fig.add_subplot(111)
        self.axes.set_facecolor('black')
        
        pylab.setp(self.axes.get_xticklabels(), fontsize=8)
        pylab.setp(self.axes.get_yticklabels(), fontsize=8)
        self.canvas = FigCanvas(self, wx.ID_ANY, self.fig)
        self.zoom_range = self.axes.axvspan(0, 0, facecolor='w', alpha=0.5)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, wx.EXPAND)

        self.SetSizer(sizer)

        self.canvas.mpl_connect('button_press_event', self.on_click)

    def clear_everything(self):
        for i in xrange(self.axes.lines.__len__()):
            del self.axes.lines[0]

    def on_click(self, event):
        if event.xdata is not None and self.axes.lines.__len__() > 0:
            self.set_zoom_box(event.xdata)
            self.plot_main(event.xdata)

    def plot_init(self, channels):
        self.channel_numbers = []
        
        for channel in channels:
            colours = self.parent.control_panel.colour_pickers[channel.channel_number].GetColour()
            col = [c / 255. for c in colours]
            self.axes.plot(channel.data[0], channel.data[1], color=col)
            self.channel_numbers.append(channel.channel_number)

        self.axes.set_xlim(channels[0].data[0][0] - PLOT_OFFSET, channels[0].data[0][-1] + PLOT_OFFSET)
        y_max = max([max(c.data[1]) for c in channels])
        y_min = min(channels[0].data[1])
        y_max_scale = Y_MAX_SCALLING * y_max
        y_min_scale = y_min - (y_max_scale - y_max)
        if y_min_scale < 0:
            y_min_scale = 0
        self.axes.set_ylim(y_min_scale, y_max_scale)
        left = channels[0].data[0][0] - DEFAULT_PLOT_RANGE / 2
        right = channels[0].data[0][0] + DEFAULT_PLOT_RANGE / 2
        self.zoom_range.set_xy([[left, 0], [right, 0], [right, 1], [left, 1]])
        self.canvas.draw()
        self.plot_main()

    def plot_main(self, start=None):
        self.parent.plot(start)

    def set_zoom_box(self, x):
        left = x - DEFAULT_PLOT_RANGE / 2
        right = x + DEFAULT_PLOT_RANGE / 2
        self.zoom_range.set_xy([[left, 0], [right, 0], [right, 1], [left, 1]])
        self.canvas.draw()

    def toggle_line_colour(self, channel, colour):
        col = [c / 255. for c in colour]
        self.axes.lines[self.channel_numbers.index(channel)].set_color(col)
        self.canvas.draw()

class MainFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, 'OPTIMA', size=(1200, 700))

        self.data_frame = DataFrame()
        self.lines = {}

        self.dpi = 100
        self.fig = Figure((1.0, 1.0), dpi=self.dpi)

        self.fig.subplots_adjust(left=0.05, right=0.97, bottom=0.05, top=0.95)

        self.axes = self.fig.add_subplot(111)
        self.axes.set_facecolor('black')
        
        pylab.setp(self.axes.get_xticklabels(), fontsize=8)
        pylab.setp(self.axes.get_yticklabels(), fontsize=8)
        self.canvas = FigCanvas(self, wx.ID_ANY, self.fig)

        self.dark_lines = []
        for i in xrange(4):
            dark_line = self.axes.axvspan(0, 0, color='r')
            self.dark_lines.append(dark_line)
            dark_line.set_visible(False)

        self.flat_lines = []
        for i in xrange(4):
            flat_line = self.axes.axvspan(0, 0, color='g')
            self.flat_lines.append(flat_line)
            flat_line.set_visible(False)

        sizer = wx.BoxSizer(wx.VERTICAL)
        self.mini_map = MiniMapPanel(self)
        sizer.Add(self.mini_map, 0, wx.GROW)
        sizer.Add(self.canvas, 1, wx.GROW)
        self.control_panel = ControlPanel(self)
        sizer.Add(self.control_panel, 0, wx.GROW)

        self.SetSizer(sizer)

    def clear_everything(self):
        for i in xrange(self.axes.lines.__len__()):
            del self.axes.lines[0]

    def open_files(self, file_list):
        self.data_frame.clear_everything()
        del self.data_frame
        self.mini_map.clear_everything()
        self.data_frame = DataFrame()
        self.data_frame.open_files(file_list)
        self.control_panel.enable(self.data_frame.channels.keys())
        self.mini_map.plot_init(self.data_frame.channels.values())

    def plot(self, plot_start=None):
        self.clear_everything()

        mi, ma = self.data_frame.get_range()

        if plot_start is None:
            plot_start = mi

          
        plot_start -= DEFAULT_PLOT_RANGE / 2 
        plot_end = plot_start + DEFAULT_PLOT_RANGE

        data_start = int(plot_start - mi)
        if data_start < 0:
            data_start = 0
        data_end = int(plot_end - mi)

        channels = self.data_frame.channels.values()

        visible_lines = [c for c in channels if self.control_panel.checkboxes[c.channel_number].GetValue()] or channels

        y_min = min([min(c.data[1][data_start:data_end + 1]) for c in visible_lines])
        y_max = max([max(c.data[1][data_start:data_end + 1]) for c in visible_lines])

        y_max += 0.1 * (y_max - y_min)
        y_min -= 0.1 * (y_max - y_min)

        self.channel_numbers = []

        for channel in channels:
            colours = self.control_panel.colour_pickers[channel.channel_number].GetColour()
            col = [c / 255. for c in colours]
            self.axes.plot(channel.data[0][data_start:data_end], channel.data[1][data_start:data_end], color=col)
            self.channel_numbers.append(channel.channel_number)
            visible = self.control_panel.checkboxes[channel.channel_number].GetValue()
            self.axes.lines[-1].set_visible(visible)

        plot_x_min = channels[0].data[0][0]
        plot_x_max = channels[0].data[0][-1]

        for i, d in enumerate(self.dark_lines):
            if not d.get_visible():
                d.set_visible(True)
            x = self.control_panel.darks[i].GetValue()
            d.set_xy([[x, 0], [x, 1]])

        for i, f in enumerate(self.flat_lines):
            if not f.get_visible():
                f.set_visible(True)
            x = self.control_panel.flats[i].GetValue()
            f.set_xy([[x, 0], [x, 1]])

        self.axes.set_xlim(plot_start, plot_end)
        self.axes.set_ylim(y_min, y_max)
        self.canvas.draw()

    def set_dark_line(self, n, v):
        x_min, x_max = self.axes.get_xlim()
        if v < x_min + DARK_LINE_OFFSET:
            plot_start = v - DARK_LINE_OFFSET + DEFAULT_PLOT_RANGE / 2
            self.plot(plot_start=plot_start)
            self.mini_map.set_zoom_box(plot_start)
        elif v > x_max - DARK_LINE_OFFSET:
            plot_start = v + DARK_LINE_OFFSET - DEFAULT_PLOT_RANGE / 2
            self.plot(plot_start=plot_start)
            self.mini_map.set_zoom_box(plot_start)
        d = self.dark_lines[n]
        d.set_xy([[v, 0], [v, 1]])
        self.canvas.draw()

    def set_flat_line(self, n, v):
        x_min, x_max = self.axes.get_xlim()
        if v < x_min + FLAT_LINE_OFFSET:
            plot_start = v - FLAT_LINE_OFFSET + DEFAULT_PLOT_RANGE / 2
            self.plot(plot_start=plot_start)
            self.mini_map.set_zoom_box(plot_start)
        elif v > x_max - FLAT_LINE_OFFSET:
            plot_start = v + FLAT_LINE_OFFSET - DEFAULT_PLOT_RANGE / 2
            self.plot(plot_start=plot_start)
            self.mini_map.set_zoom_box(plot_start)
        f = self.flat_lines[n]
        f.set_xy([[v, 0], [v, 1]])
        self.canvas.draw()

    def toggle_line(self, channel, visible=True):
        self.axes.lines[self.channel_numbers.index(channel)].set_visible(visible)
        #self.canvas.draw()
        x_min, x_max = self.axes.get_xlim()
        plot_start = (x_max + x_min) / 2
        self.plot(plot_start=plot_start)

    def toggle_line_colour(self, channel, colour):
        col = [c / 255. for c in colour]
        self.axes.lines[self.channel_numbers.index(channel)].set_color(col)
        self.canvas.draw()

class OpenFilesDialog(wx.FileDialog):
    def __init__(self, default_dir):
        message = 'open files'
        wx.FileDialog.__init__(self, parent=None, message=message, defaultDir=default_dir, wildcard=self.supported_files(), style=wx.OPEN|wx.FD_MULTIPLE)

    def supported_files(self):
        return 'FITS file (*.fits;*.FITS)|*.fits;*.FITS'

class PlotLC(wx.Dialog):
    def __init__(self, lc):
        wx.Dialog.__init__(self, parent=None, id=-1, title='Lightcurve plot')
        sizer = wx.BoxSizer(wx.VERTICAL)

        self.dpi = 100
        self.fig = Figure(dpi=self.dpi)
        self.canvas = FigCanvas(self, -1, self.fig)

        self.axes = self.fig.add_subplot(111)
        self.toolbar = NavigationToolbar(self.canvas)

        self.axes.plot(lc)
        self.canvas.draw()

        button = wx.Button(self, wx.ID_OK)

        sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        sizer.Add(self.toolbar, 0, wx.EXPAND)
        sizer.Add(button)

        self.SetSizer(sizer)
        sizer.Fit(self)

class SaveFileDialog(wx.FileDialog):
    def __init__(self, default_dir, default_file, data_type, extension):
        message = 'save %s %s file' % (data_type, extension)
        default_file += '_' + data_type + '.' + extension
        wx.FileDialog.__init__(self, parent=None, message=message, defaultDir=default_dir, defaultFile=default_file, wildcard=self.supported_files(extension), style=wx.SAVE|wx.OVERWRITE_PROMPT)

    def supported_files(self, extension):
        EXTENSIONS = {'fits': 'FITS file (*.fits;*.FITS)|*.fits;*.FITS',
                        'txt': 'text file (*.txt)|*.txt',
                    }
        return EXTENSIONS[extension]

if __name__ == '__main__':
    app = wx.App()
    app.frame = MainFrame()
    app.frame.Show()
    app.MainLoop()
