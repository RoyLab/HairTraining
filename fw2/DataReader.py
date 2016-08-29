from common_tools import *
from coordinates import *
import array
import numpy as np

class FrameData:
    def __init__(self):
        self.headMotion = None
        self.position = None
        self.direction = None

class HairHeader:
    def __init__(self):
        self.nParticle = None
        self.nHair = None
        self.factor = None

class HairDataReader:
    def __init__(self, fileName, args):
        self.type = tryGetPara(args, "type")
        if self.type == "anim2":
            self.file = open(fileName, 'rb')
            self.nFrame = readInt(self.file)
            self.nParticle = readInt(self.file)
            self.start = self.file.tell()
            self.offset = 17*4 + 6*self.nParticle
            self.frameId = None
            self.rewind()

    def rewind(self):
        if self.type == "anim2":
            self.file.seek(self.start)
        self.frameId = 0

    def curFrame(self):
        return self.frameId

    def getNextFrame(self):
        self.frameId += 1
        if self.frameId >= self.nFrame:
            self.rewind()

        if self.type == "anim2":
            readInt(self.file)
            frame = FrameData()

            tmp = array.array('f')
            tmp.fromfile(self.file, 16)
            frame.headMotion = MatrixToRt(np.matrix(tmp).reshape((4,4)))

            frame.position = array.array('f')
            frame.position.fromfile(self.file, self.nParticle*3)

            frame.direction = array.array('f')
            frame.direction.fromfile(self.file, self.nParticle*3)

            return frame

    def randomGetFrame(self, n):
        if self.type == "anim2":
            tmp = self.start + self.offset * n
            self.file.seek(tmp)

        self.frameId = n-1
        self.getNextFrame()

if __name__== "__main__":
    fileName = r"D:\Data\c0524\c0514.anim2"
    a = HairDataReader(fileName)