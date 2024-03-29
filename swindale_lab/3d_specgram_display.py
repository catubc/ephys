"""PyQt OpenGL example, modified from PySide example, which was a port of the opengl/hellogl
example from Qt v4.x"""

from __future__ import division

import sys
import time

# instantiate an IPython embedded shell which shows up in the terminal on demand
# and on every exception:
from IPython.terminal.ipapp import load_default_config
from IPython.terminal.embed import InteractiveShellEmbed
config = load_default_config()
# automatically call the pdb debugger after every exception, override default config:
config.TerminalInteractiveShell.pdb = True
ipshell = InteractiveShellEmbed(display_banner=False, config=config)

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt

from OpenGL import GL, GLU
import numpy as np

RED = 255, 0, 0
ORANGE = 255, 127, 0
YELLOW = 255, 255, 0
GREEN = 0, 255, 0
CYAN = 0, 255, 255
LIGHTBLUE = 0, 127, 255
BLUE = 0, 0, 255
VIOLET = 127, 0, 255
MAGENTA = 255, 0, 255
GREY = 85, 85, 85
WHITE = 255, 255, 255
CMAP = np.array([RED, ORANGE, YELLOW, GREEN, CYAN, LIGHTBLUE, BLUE, VIOLET, MAGENTA,
                 GREY, WHITE], dtype=np.uint8)

'''
def normdeg(angle):
    return angle % 360
'''

class Window(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.glWidget = GLWidget(parent=self)
        mainLayout = QtGui.QHBoxLayout()
        mainLayout.addWidget(self.glWidget)
        self.setLayout(mainLayout)
        self.setWindowTitle(self.tr("OpenGL test"))

    def keyPressEvent(self, event):
        self.glWidget.keyPressEvent(event) # pass it down


class GLWidget(QtOpenGL.QGLWidget):
    def __init__(self, parent=None):
        QtOpenGL.QGLWidget.__init__(self, parent)
        self.lastPos = QtCore.QPoint()
        self.focus = np.float32([0, 0, 0]) # init camera focus

        format = QtOpenGL.QGLFormat()
        #format.setVersion(3, 0) # not available in PyQt 4.7.4
        # set to color index mode, unsupported in OpenGL >= 3.1, don't know how to load
        # GL_ARB_compatibility extension, and for now, can't force OpenGL 3.0 mode.
        # Gives "QGLContext::makeCurrent(): Cannot make invalid context current." error:
        #format.setRgba(False)
        format.setDoubleBuffer(True) # works fine
        self.setFormat(format)
        #QtOpenGL.QGLFormat.setDefaultFormat(format)
        '''
        c = QtGui.qRgb
        cmap = [c(255, 0, 0), c(0, 255, 0), c(0, 0, 255), c(255, 255, 0), c(255, 0, 255)]
        colormap = QtOpenGL.QGLColormap()
        colormap.setEntries(cmap)
        self.setColormap(colormap)
        '''
        
        #Pick random points
        if False:
            self.npoints = 100000
            if self.npoints > 2**24-2: # the last one is the full white bg used as a no hit
                raise OverflowError("Can't pick from more than 2**24-2 sids")
            self.points = np.float32(np.random.random((self.npoints, 3))) - 0.5
            self.sids = np.arange(self.npoints)
            self.nids = self.sids % len(CMAP)
            self.colors = CMAP[self.nids] # uint8

            alpha_array = np.zeros(len(self.colors),dtype=np.int16)
            self.colors = np.insert(np.float32(self.colors), 3,values=0.5, axis=1)
        
        else:
            
            #dir_name = '/media/cat/12TB/in_vivo/tim/dongsheng/2015-11-27/2015-11-27-2-10electrodein-iso1/2015-11-27-2-10electrodein-iso1_raw.tsf'
            #dir_name = '/media/cat/12TB/in_vivo/nick/ptc21/61-tr5c-blankscreen/61-tr5c-blankscreen.lfp.zip'
            dir_name = '/media/cat/12TB/in_vivo/nick/ptc17/44-tr2b-blankscreen/44-tr2b-blankscreen.lfp.zip'
            self.points = []
            self.colors = []
            for h in range(10):
                data = np.load(dir_name+"_envelope_"+str((h+1)*10)+"hz.npy")

                for k in range(len(data)):
                    print "Loading ch: ", k
                    npoints = len(data[k])

                    #Load origin of points
                    time_points = np.zeros((npoints, 3), dtype=np.float32) 
                    time_points[:,2]=np.arange(npoints)*2    #Insert depth direction in the z axis (3rd column)
                    
                    #Compute Amplitudes
                    amplitudes = np.float32(data[k])*1E-3
                    amp0 = np.zeros((npoints,3), dtype=np.float32); amp1 = amp0.copy(); amp2 = amp0.copy(); amp3=amp0.copy(); amp4 = amp0.copy()
                    amp1[:,0] = 0+h*20; amp1[:,1] = amplitudes*0.5+5*k; amp1[:,2] = amplitudes*0.5
                    amp2[:,0] = 0+h*20; amp2[:,1] = amplitudes*0.5+5*k; amp2[:,2] = amplitudes*(-0.5)
                    amp3[:,0] = 0+h*20; amp3[:,1] = amplitudes*(-0.5)+5*k; amp3[:,2] = amplitudes*(-0.5)
                    amp4[:,0] = 0+h*20; amp4[:,1] = amplitudes*(-0.5)+5*k; amp4[:,2] = amplitudes*0.5
                    
                    #Convert points to quads
                    points = np.zeros((npoints*4,3),dtype=np.float32)
                    points[::4] = time_points+amp1
                    points[1::4] = time_points+amp2
                    points[2::4] = time_points+amp3
                    points[3::4] = time_points+amp4
                    self.points.extend(points)

                    #Set colors
                    colors = [[255,255,255, 1.0]]*npoints*4 # uint8
                    self.colors.extend(colors)

            self.npoints = len(self.points)
            self.sids = np.arange(self.npoints)

            print self.npoints
            print len(self.colors)
            self.points = np.vstack((self.points))
            print self.points[0:10]
            #quit()
            
            #alpha_array = np.zeros(len(self.colors),dtype=np.int16)
            #self.colors = np.insert(np.float32(self.colors), 3,values=0.5, axis=1)

        # encode sids in RGB - used for picking points 
        r = self.sids // 256**2
        rem = self.sids % 256**2 # remainder
        g = rem // 256
        b = rem % 256
        self.rgbsids = np.zeros((self.npoints, 3), dtype=np.uint8)
        self.rgbsids[:, 0] = r
        self.rgbsids[:, 1] = g
        self.rgbsids[:, 2] = b
        #print self.rgbsids
        #print self.rgbsids.dtype

    def minimumSizeHint(self):
        return QtCore.QSize(50, 50)

    def sizeHint(self):
        return QtCore.QSize(600, 400)

    def initializeGL(self):
        GL.glClearColor(0.0, 0.0, 0.0, 1.0) # same as default
        GL.glClearDepth(1.0) # same as default
        GL.glEnable(GL.GL_DEPTH_TEST) # display points according to occlusion, not order of plotting
        #GL.glEnable(GL.GL_POINT_SMOOTH) # doesn't seem to work right, proper way to antialiase?
        #GL.glEnable(GL.GL_LINE_SMOOTH) # works better
        #GL.glPointSize(1.5) # truncs to the nearest pixel if antialiasing is off
        #GL.glShadeModel(GL.GL_FLAT)
        #GL.glEnable(GL.GL_CULL_FACE) # only useful for solids
        GL.glTranslate(0, 0, -3) # init camera distance from origin

        #Enable alpha transparency
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)

    def paintGL(self):
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

        # Don't load identity matrix. Do all transforms in place against current matrix
        # and take advantage of OpenGL's state-machineness.
        # Sure, you might get round-off error over time, but who cares? If you wanna return
        # to a specific focal point on 'f', that's when you really need to first load the
        # identity matrix before doing the transforms
        #GL.glLoadIdentity() # loads identity matrix into top of matrix stack

        # viewing transform for camera: where placed, where pointed, which way is up:
        #GLU.gluLookAt()
        #GL.glScale() # modelling transformation, lets you stretch your objects
        
        GL.glEnableClientState(GL.GL_COLOR_ARRAY);
        GL.glEnableClientState(GL.GL_VERTEX_ARRAY);

        #******* Draw points *******
        #sizes = 5 #np.zeros(self.npoints,dtype=np.int16)+10
        #GL.glPointSize(sizes)             #;//set point size to 10 pixels
        ##GL.glColorPointerub(self.colors) # unsigned byte, ie uint8
        #GL.glColorPointer(4, GL.GL_FLOAT, 0, self.colors)       #colors contain alpha transparency values also now
        #GL.glVertexPointerf(self.points) # float32
        #GL.glDrawArrays(GL.GL_POINTS, 0, self.npoints)

        #****** Draw quads ******
        GL.glColorPointer(4, GL.GL_FLOAT, 0, self.colors) # unsigned byte, ie uint8
        GL.glVertexPointerf(self.points) # float32
        GL.glDrawArrays(GL.GL_QUADS, 0, len(self.points))
        GL.glColor3ub(255, 255, 255)

        # might consider using buffer objects for even more speed (less unnecessary vertex
        # data from ram to vram, I think). Apparently, buffer objects don't work with
        # color arrays?

        #GL.glFlush() # forces drawing to begin, only makes difference for client-server?
        self.swapBuffers() # doesn't seem to be necessary, even though I'm in double-buffered
                           # mode with the back buffer for RGB sid encoding, but do it anyway
                           # for completeness

        # print the modelview matrix
        #print(self.MV)

    def resizeGL(self, width, height):
        GL.glViewport(0, 0, width, height)
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()
        # fov (deg) controls amount of perspective, and as a side effect initial apparent size
        GLU.gluPerspective(45, width/height, 0.0001, 1000) # fov, aspect, nearz & farz
                                                           # clip planes
        GL.glMatrixMode(GL.GL_MODELVIEW)

    def get_MV(self):
        """Return modelview matrix"""
        return GL.glGetDoublev(GL.GL_MODELVIEW_MATRIX) # I think this acts like a copy

    def set_MV(self, MV):
        GL.glLoadMatrixd(MV)

    MV = property(get_MV, set_MV)

    # modelview matrix is column major, so we work on columns instead of rows
    def getViewRight(self):
        """View right vector: 1st col of modelview matrix"""
        return self.MV[:3, 0]

    def getViewUp(self):
        """View up vector: 2nd col of modelview matrix"""
        return self.MV[:3, 1]

    def getViewNormal(self):
        """View normal vector: 3rd col of modelview matrix"""
        return self.MV[:3, 2]

    def getTranslation(self):
        """Translation vector: 4th row of modelview matrix"""
        return self.MV[3, :3]

    def setTranslation(self, vt):
        """Translation vector: 4th row of modelview matrix"""
        MV = self.MV
        MV[3, :3] = vt
        self.MV = MV

    def getDistance(self):
        v = self.getTranslation()
        #return np.sqrt((v**2).sum()) # from data origin
        return np.sqrt(((v-self.focus)**2).sum()) # from focus

    def pan(self, dx, dy):
        """Translate along view right and view up vectors"""
        d = self.getDistance()
        vr = self.getViewRight()
        vr *= dx*d
        GL.glTranslate(vr[0], vr[1], vr[2])
        vu = self.getViewUp()
        vu *= dy*d
        GL.glTranslate(vu[0], vu[1], vu[2])

    def zoom(self, dr):
        """Translate along view normal vector"""
        d = self.getDistance()
        vn = self.getViewNormal()
        vn *= dr*d
        GL.glTranslate(vn[0], vn[1], vn[2])

    def pitch(self, dangle): # aka elevation
        """Rotate around view right vector"""
        vr = self.getViewRight()
        GL.glTranslate(*self.focus)
        GL.glRotate(dangle, *vr)
        GL.glTranslate(*-self.focus)

    def yaw(self, dangle): # aka azimuth
        """Rotate around view up vector"""
        vu = self.getViewUp()
        GL.glTranslate(*self.focus)
        GL.glRotate(dangle, *vu)
        GL.glTranslate(*-self.focus)

    def roll(self, dangle):
        """Rotate around view normal vector"""
        vn = self.getViewNormal()
        GL.glTranslate(*self.focus)
        GL.glRotate(dangle, *vn)
        GL.glTranslate(*-self.focus)

    def panTo(self, p=None):
        """Translate along view right and view up vectors such that data point p is
        centered in the viewport. Not entirely sure why or how this works, figured
        it out using guess and test"""
        if p == None:
            p = self.focus
        MV = self.MV
        vr = self.getViewRight()
        vu = self.getViewUp()
        p = -p
        x = np.dot(p, vr) # dot product
        y = np.dot(p, vu)
        MV[3, :2] = x, y # set first two entries of 4th row to x, y
        self.MV = MV

    def pick(self):
        globalPos = QtGui.QCursor.pos()
        pos = self.mapFromGlobal(globalPos)
        width = self.size().width()
        height = self.size().height()
        x = pos.x()
        y = height - pos.y()
        if not (0 <= x < width and 0 <= y < height):
            print('cursor out of range')
            return
        '''
        # for speed, 1st check if there are any non-black pixels around cursor:
        GL.glReadBuffer(GL.GL_FRONT)
        frontbuffer = GL.glReadPixelsub(0, 0, width, height, GL.GL_RGB) # unsigned byte
        #rgb = frontbuffer[y-1:y+2, x-1:x+2] # +/- 1 pix
        rgb = frontbuffer[y, x]
        #print('frontbuffer:')
        #print rgb
        if (rgb == 0).all():
            print('nothing to return')
            return # nothing to pick
        '''
        # drawing encoded RGB values to back buffer
        #GL.glDrawBuffer(GL_BACK) # shouldn't be necessary, defaults to back
        GL.glClearColor(1.0, 1.0, 1.0, 1.0) # highest possible RGB means no hit
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

        GL.glEnableClientState(GL.GL_COLOR_ARRAY);
        GL.glEnableClientState(GL.GL_VERTEX_ARRAY);
        GL.glColorPointerub(self.rgbsids) # unsigned byte, ie uint8
        GL.glVertexPointerf(self.points) # float32

        GL.glDrawArrays(GL.GL_POINTS, 0, self.npoints) # to back buffer

        GL.glClearColor(0.0, 0.0, 0.0, 1.0) # restore to default black

        # grab back buffer
        GL.glReadBuffer(GL.GL_BACK)

        # find rgb at cursor coords, decode sid
        # unsigned byte, x, y is bottom left:
        backbuffer = GL.glReadPixelsub(x, y, 1, 1, GL.GL_RGB)
        rgb = backbuffer[0, 0]
        r, g, b = rgb
        sid = r*256**2 + g*256 + b
        if sid == 2**24 - 1:
            #print('no hit')
            return
        nid = sid % len(CMAP)
        color = CMAP[nid]
        print('backbuffer.shape: %r' % (backbuffer.shape,))
        print rgb
        print('sid, nid, color: %d, %d, %r' % (sid, nid, color))
        '''
        # TODO: this isn't even necessary, since back buffer gets overdrawn anyway on next
        # updateGL():
        # restore front buffer to back
        GL.glReadBuffer(GL.GL_FRONT)
        #GL.glDrawBuffer(GL_BACK) # shouldn't be necessary, defaults to back
        GL.glRasterPos(0, 0) # destination position in draw buffer
        GL.glCopyPixels(0, 0, width, height, GL.GL_COLOR) # copies from read to draw buffer
        '''
        #self.swapBuffers() # don't even need to swap buffers, cuz we haven't changed the scene
        return sid

    def mousePressEvent(self, event):
        self.lastPos = QtCore.QPoint(event.pos())

    def mouseMoveEvent(self, event):
        buttons = event.buttons()
        modifiers = event.modifiers()
        dx = event.x() - self.lastPos.x()
        dy = event.y() - self.lastPos.y()

        if buttons == QtCore.Qt.LeftButton:
            if modifiers == Qt.ControlModifier:
                self.roll(-0.5*dx - 0.5*dy)
            elif modifiers == Qt.ShiftModifier:
                self.pan(dx/600, -dy/600) # qt viewport y axis points down
            else:
                self.yaw(0.5*dx)
                self.pitch(0.5*dy)
        elif buttons == QtCore.Qt.RightButton:
            self.zoom(-dy/500) # qt viewport y axis points down

        self.updateGL()
        self.lastPos = QtCore.QPoint(event.pos())

    def wheelEvent(self, event):
        self.zoom(event.delta() / 1000)
        self.updateGL()

    def keyPressEvent(self, event):
        key = event.key()
        modifiers = event.modifiers()

        if modifiers == Qt.ControlModifier:
            if key == Qt.Key_Left:
                self.roll(5)
            elif key == Qt.Key_Right:
                self.roll(-5)
            elif key == Qt.Key_Up:
                self.zoom(0.05)
            elif key == Qt.Key_Down:
                self.zoom(-0.05)
        elif modifiers == Qt.ShiftModifier:
            if key == Qt.Key_Left:
                self.pan(-0.05, 0)
            elif key == Qt.Key_Right:
                self.pan(0.05, 0)
            elif key == Qt.Key_Up:
                self.pan(0, 0.05)
            elif key == Qt.Key_Down:
                self.pan(0, -0.05)
        else:
            if key == Qt.Key_Left:
                self.yaw(-5)
            elif key == Qt.Key_Right:
                self.yaw(5)
            elif key == Qt.Key_Up:
                self.pitch(-5)
            elif key == Qt.Key_Down:
                self.pitch(5)
            elif key == Qt.Key_0: # reset focus to origin
                self.focus = np.float32([0, 0, 0])
                self.panTo() # pan to new focus
            elif key == Qt.Key_F: # reset focus to cursor position
                sid = self.pick()
                if sid != None:
                    self.focus = self.points[sid]
                    self.panTo() # pan to new focus
            elif key == Qt.Key_P:
                self.pick()
            elif key == Qt.Key_S:
                ipshell()

        self.updateGL()

    '''
    # this specifies a display list, which is sent once, compiled, and then simply referenced
    # later every time the display needs to be updated. However, display lists are static once
    # compiled - none of their attributes can be changed
    def makeObject(self):
        genList = GL.glGenLists(1)
        GL.glNewList(genList, GL.GL_COMPILE)

        GL.glBegin(GL.GL_QUADS)

        x1 = +0.06
        y1 = -0.14
        x2 = +0.14
        y2 = -0.06
        x3 = +0.08
        y3 = +0.00
        x4 = +0.30
        y4 = +0.22

        self.quad(x1, y1, x2, y2, y2, x2, y1, x1)
        self.quad(x3, y3, x4, y4, y4, x4, y3, x3)

        self.extrude(x1, y1, x2, y2)
        self.extrude(x2, y2, y2, x2)
        self.extrude(y2, x2, y1, x1)
        self.extrude(y1, x1, x1, y1)
        self.extrude(x3, y3, x4, y4)
        self.extrude(x4, y4, y4, x4)
        self.extrude(y4, x4, y3, x3)

        Pi = 3.14159265358979323846
        NumSectors = 200

        for i in range(NumSectors):
            angle1 = (i * 2 * Pi) / NumSectors
            x5 = 0.30 * math.sin(angle1)
            y5 = 0.30 * math.cos(angle1)
            x6 = 0.20 * math.sin(angle1)
            y6 = 0.20 * math.cos(angle1)

            angle2 = ((i + 1) * 2 * Pi) / NumSectors
            x7 = 0.20 * math.sin(angle2)
            y7 = 0.20 * math.cos(angle2)
            x8 = 0.30 * math.sin(angle2)
            y8 = 0.30 * math.cos(angle2)

            self.quad(x5, y5, x6, y6, x7, y7, x8, y8)

            self.extrude(x6, y6, x7, y7)
            self.extrude(x8, y8, x5, y5)

        GL.glEnd()
        GL.glEndList()

        return genList

    def quad(self, x1, y1, x2, y2, x3, y3, x4, y4):
        self.qglColor(self.trolltechGreen)

        GL.glVertex3d(x1, y1, -0.05)
        GL.glVertex3d(x2, y2, -0.05)
        GL.glVertex3d(x3, y3, -0.05)
        GL.glVertex3d(x4, y4, -0.05)

        GL.glVertex3d(x4, y4, +0.05)
        GL.glVertex3d(x3, y3, +0.05)
        GL.glVertex3d(x2, y2, +0.05)
        GL.glVertex3d(x1, y1, +0.05)

    def extrude(self, x1, y1, x2, y2):
        self.qglColor(self.trolltechGreen.darker(250 + int(100 * x1)))

        GL.glVertex3d(x1, y1, +0.05)
        GL.glVertex3d(x2, y2, +0.05)
        GL.glVertex3d(x2, y2, -0.05)
        GL.glVertex3d(x1, y1, -0.05)
    '''
if __name__ == '__main__':
    # prevents "The event loop is already running" messages when calling ipshell():
    QtCore.pyqtRemoveInputHook()
    app = QtGui.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())

