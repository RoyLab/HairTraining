from frame import Frame
from progressbar import *
from GraphBuilder import *
import numpy as np
import time

from multiprocessing import Pool, TimeoutError, Queue, Manager

MAX_PROCESSERS = 3

class ParallelHooker(object):
    def __init__(self, number, post):
        self.nFrame = number
        self.i = -1
        self.manager = Manager()
        self.queue = self.manager.Queue(1)
        self.output = self.manager.Queue()
        self.pool = Pool(processes=MAX_PROCESSERS)

        self.workers = MAX_PROCESSERS
        self.nInputed = 0
        self.nAdded = 0
        self.nOutputed = 0

        self.postFunc = post

        return

    def startLoop(self, title="start loop:"):
        print title
        self.bar = ProgressBar().start()
        return

    def endLoop(self):
        while self.nAdded < self.nFrame:
            self.pool.apply_async(self.handle, ())
            self.nAdded += 1

        self.pool.close()
        while self.nOutputed < self.nFrame:
            res = self.output.get()
            self.handle(res)

            self.i += 1
            self.bar.update((self.i+1)*100/self.nFrame)
            self.nOutputed += 1

        self.pool.join()
        self.bar.finish()

    def newFrame(self):
        self.frame = Frame()
        return

    def postFrame(self):

        # check if all workers is used, else block, wait workers
        while self.workers < 1:
            while not self.output.empty() and self.nOutputed < self.nFrame:
                res = self.output.get()
                self.handle(res)

                self.i += 1
                self.bar.update((self.i+1)*100/self.nFrame)
                self.nOutputed += 1
                self.workers += 1

            # keep waiting workers be freed
            if self.workers < 1:
                time.sleep(1)

        # now there are workers
        # empty the queue by adding tasks
        # if no workers left, then leave the data and go on
        # it will be blocked when putting data
        while self.workers > 0 and not self.queue.empty()\
         and self.nAdded < self.nFrame:
            self.pool.apply_async(self.handle, (self.queue, self.output))
            self.nAdded += 1
            self.workers -= 1

        self.queue.put(self.frame, True)

        return

    def handle(self, res):
        return

    def dataHooker(self, name, sz, arr):
        self.frame.loadIntoMemory(name, sz, arr)
        return
