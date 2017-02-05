
"""A skeletal copy of classes from dimstim required for interpreting
the textheader generated by dimstim.__version__ >= 0.16"""

from __future__ import division
from __future__ import print_function

import os
import struct
from copy import copy
import cStringIO
from struct import Struct, unpack

import numpy as np

import core
from core import dictattr, TAB

class Experiment(object):
    """Base Experiment class, all dimstim experiments inherit from this"""
    def __init__(self, script=None, static=None, dynamic=None, variables=None, runs=None,
                 blanksweeps=None):
        self.script = script # Experiment script file name
        self.static = static # StaticParams object
        self.dynamic = dynamic # DynamicParams object
        self.variables = variables # Variables object
        self.runs = runs # Runs object
        self.blanksweeps = blanksweeps # BlankSweeps object

class Movie(Experiment):
    
    def load(self, file_name, asarray=False, flip=False):
        """Load movie frames"""
        # figure out the local path to the same movie file:
        #pathparts = core.pathdecomp(self.static.fname) # as it existed on the stim computer
        #movi = pathparts.index('mov')
        #tail = os.path.join(pathparts[movi+1:]) # everything after 'mov' folder
        #MOVIEPATH = get_ipython().user_ns['MOVIEPATH']
        #fullfname = os.path.join(MOVIEPATH, *tail) # full fname with local MOVIEPATH
        self.f = file(file_name, 'rb') # open the movie file for reading in binary format
        headerstring = self.f.read(5)
        if headerstring == 'movie': # a header has been added to the start of the file
            self.ncellswide, = struct.unpack('H', self.f.read(2)) # 'H'== unsigned short int
            self.ncellshigh, = struct.unpack('H', self.f.read(2))
            self.nframes, = struct.unpack('H', self.f.read(2))
            if self.nframes == 0:
                # this was used in ptc15 mseq movies to indicate 2**16 frames, shouldn't
                # really worry about this, cuz we're using slightly modified mseq movies now
                # that don't have the extra frame at the end that the ptc15 movies had (see
                # comment in Experiment module), and therefore never have a need to indicate
                # 2**16 frames
                self.nframes = 2**16
            self.offset = self.f.tell() # header is 11 bytes long
        else: # there's no header at the start of the file, set the file pointer back to the beginning and use these hard coded values:
            self.f.seek(0)
            self.ncellswide = self.ncellshigh = 64
            self.nframes = 6000
            self.offset = self.f.tell() # header is 0 bytes long
        self.framesize = self.ncellshigh*self.ncellswide

        # read in all of the frames. Maybe check first to see if file is > 1GB. If so,
        # _loadaslist() to prevent trying to allocate one huge piece of contiguous memory and
        # raising a MemoryError, or worse, segfaulting
        if asarray:
            self._loadasarray(flip)
        else:
            self._loadaslist(flip)
        leftover = self.f.read() # check if there are any leftover bytes in the file
        if leftover != '':
            pprint(leftover)
            print(self.ncellswide, self.ncellshigh, self.nframes)
            raise RuntimeError('There are unread bytes in movie file %r. Width, height, or nframes is incorrect in the movie file header.' % self.static.fname)
        self.f.close() # close the movie file

    def _loadasarray(self, flip):
        self.frames = np.fromfile(self.f, dtype=np.uint8, count=self.nframes*self.framesize)
        self.frames.shape = (self.nframes, self.ncellshigh, self.ncellswide)
        if flip:
            self.frames = self.frames[::, ::-1, ::] # flip all frames vertically for OpenGL's bottom left origin

    def _loadaslist(self, flip):
        self.frames = []
        for framei in xrange(self.nframes): # one frame at a time...
            frame = np.fromfile(self.f, dtype=np.uint8, count=self.framesize) # load the next frame
            frame.shape = (self.ncellshigh, self.ncellswide)
            if flip:
                frame = frame[::-1, ::] # flip all frames vertically for OpenGL's bottom left origin
            self.frames.append(frame)



class StimulusHeader(object):
    """Stimulus display header"""
    # Stimulus header constants
    OLD_STIMULUS_HEADER_FILENAME_LEN = 16
    STIMULUS_HEADER_FILENAME_LEN = 64
    NVS_PARAM_LEN = 749
    PYTHON_TBL_LEN = 50000

    def __init__(self):
        
        self.version = 110

    def __len__(self):
        if self.version == 100: # Cat < 15
            return 4 + self.OLD_STIMULUS_HEADER_FILENAME_LEN + self.NVS_PARAM_LEN*4 + 28
        elif self.version == 110: # Cat >= 15
            return (4 + self.STIMULUS_HEADER_FILENAME_LEN + self.NVS_PARAM_LEN*4 +
                    self.PYTHON_TBL_LEN + 28)

    def parse(self, f):
        self.header = f.read(2).rstrip(NULL) # always 'DS'?
        self.version, = unpack('H', f.read(2))
        if self.version not in (100, 110): # Cat < 15, Cat >= 15
            raise ValueError, 'Unknown stimulus header version %d' % self.version
        if self.version == 100: # Cat < 15 has filename field length == 16
            # ends with a NULL followed by spaces for some reason,
            # at least in Cat 13 file 03 - ptumap#751a_track5_m-seq.srf
            self.filename = f.read(self.OLD_STIMULUS_HEADER_FILENAME_LEN).rstrip().rstrip(NULL)
        elif self.version == 110: # Cat >= 15 has filename field length == 64
            self.filename = f.read(self.STIMULUS_HEADER_FILENAME_LEN).rstrip(NULL)
        # NVS binary header, array of single floats
        self.parameter_tbl = list(unpack('f'*self.NVS_PARAM_LEN, f.read(4*self.NVS_PARAM_LEN)))
        for parami, param in enumerate(self.parameter_tbl):
            if str(param) == '1.#QNAN':
                # replace 'Quiet NAN' floats with Nones. This won't work for Cat < 15
                # because NVS display left empty fields as NULL instead of NAN
                self.parameter_tbl[parami] = None
        # dimstim's text header
        if self.version == 110: # only Cat >= 15 has the text header
            self.python_tbl = f.read(self.PYTHON_TBL_LEN).rstrip()
        self.screen_width, = unpack('f', f.read(4)) # cm, single float
        self.screen_height, = unpack('f', f.read(4)) # cm
        self.view_distance, = unpack('f', f.read(4)) # cm
        self.frame_rate, = unpack('f', f.read(4)) # Hz
        self.gamma_correct, = unpack('f', f.read(4))
        self.gamma_offset, = unpack('f', f.read(4))
        self.est_runtime, = unpack('H', f.read(2)) # in seconds
        self.checksum, = unpack('H', f.read(2))

