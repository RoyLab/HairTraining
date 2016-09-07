import spams
import numpy as np
from DataReader import HairDataReader, HairHeader
import local_para_small as para
import coordinates as cd
import sys
from scipy.optimize import minimize
from common_tools import *
import math

import gc
import crash_on_ipy
import ipdb
import cPickle as pkl

def normalizePerColumn(X):
     X /= np.tile(np.sqrt((X*X).sum(axis=0)),(X.shape[0],1))


def evalError(x, m, lambda1):
    res = (((x * m[0]).A1.dot(x) - 2 * m[1].dot(x) + m[2])*lambda1 +
           (x * m[3]).A1.dot(x) - 2 * m[4].dot(x) + m[5])
    return res

def evalErrorGroundTruth(w, Rs, Ts, s0, gt):
    """ w: 1 x nGuide, R/Ts: nGuide * R/t """
    n = len(w)
    Rs = np.matrix(Rs.reshape(n, -1))
    Ts = np.matrix(Ts.reshape(n, -1))

    R = (w*Rs).reshape((3,3))
    T = (w*Ts).A1

    s1 = cd.point_trans(s0, R, T)
    diff = np.array(s1 - gt).flatten()
    return diff.dot(diff)


def evalDerive(x, m, lambda1):
    res = (2 * m[0].dot(x) - 2 * m[1]).A1*lambda1 + \
           (2 * m[3].dot(x) - 2 * m[4]).A1
    return res

def estimateWeight(sf, BgsfR, Bgsft, constrain):
    """sf: ndarray, nframe * ({pi}, {ti}), exclude 0, so n=factor-1
    BgsfR/t: tuple, n * nFrame * {R/t}"""
    lambda1 = 1.0

    nFrame = sf.shape[0]
    nGuide = BgsfR.shape[0]

    pmat_square = np.matrix(np.zeros((nGuide, nGuide)))
    tmat_square = np.matrix(np.zeros((nGuide, nGuide)))
    plinearmat = np.matrix(np.zeros(nGuide))
    tlinearmat = np.matrix(np.zeros(nGuide))
    pcmat_square = 0.0
    tcmat_square = 0.0


    for i in range(1, nFrame):
        pc, tc = sf[i]
        Rs = BgsfR[:, i, ...]
        Ts = Bgsft[:, i, ...]
        pos_, dir_ = cd.point_trans(sf[0], Rs, Ts, batch=True) # pos_, dir_ are flattened 1 x 3nGuide

        posMat = np.matrix(pos_.reshape(nGuide, 3))
        dirMat = np.matrix(dir_.reshape(nGuide, 3))
        pmat_square += posMat * posMat.T
        tmat_square += dirMat * dirMat.T
        plinearmat += pc * posMat.T
        tlinearmat += tc * dirMat.T
        pcmat_square += pc.dot(pc)
        tcmat_square += tc.dot(tc)

    mats = (pmat_square, plinearmat.A1, pcmat_square, tmat_square, tlinearmat.A1, tcmat_square)

    initx = np.array([1.0 / nGuide] * nGuide)
    res = minimize(evalError, initx, args=(mats, lambda1), jac=evalDerive,
                   options={'disp': True, 'ftol':1e-12, 'maxiter':100}, method='SLSQP', constraints=constrain, )

    error = evalError(res.x,mats, lambda1)

    # error = 0.0

    # for i in range(1, nFrame):
    #     Rs = BgsfR[:, i, ...]
    #     Ts = Bgsft[:, i, ...]
    #     error += evalErrorGroundTruth(res.x, Rs, Ts, sf[0], sf[i])

    return res, error

class XWrapper:
    def __init__(self, X, offset):
        self.X = X
        self.offset = offset
        self.nFrame = X.shape[0] / offset

    def getState(self, i, fId = None):
        if not fId:
            return self.X[:, i].A1
        else:
            return self.X[fId*self.offset:(fId+1)*self.offset, i].A1

def SCGetMatrixAndHeader(fileName, nFrame=None):
    reader = HairDataReader(fileName, {'type':'anim2'})
    factor = para.factor
    offset = factor * 3
    spcereg = para.lambda_balance  # regularize the space item
    if not nFrame:
        nFrame = reader.nFrame

    mat = []
    for i in xrange(nFrame):
        if i % 10 == 0:
            sys.stdout.write("\rReading frame %d..." % i)
        frame = reader.getNextFrame()
        pos = np.array(frame.position)
        dir = np.array(frame.direction)
        rigid = frame.headMotion
        invR = np.linalg.inv(rigid[0])
        pos, dir = cd.inverseRigidTrans(invR, rigid[1], pos, dir, batch=True)
        pos.shape = -1, offset
        dir.shape = -1, offset
        mat.append(pos*spcereg)
        mat.append(dir)
    print "\rFinished Reading!"

    header = HairHeader()
    header.nParticle = reader.nParticle
    header.factor = factor
    header.nHair = header.nParticle / factor

    X = np.hstack(mat).transpose()
    return X, header, XWrapper(X, offset*2)

def pickGuideHair(D, X):
    D = D.transpose()
    X = X.transpose()
    nHair = X.shape[0]
    nGuide = D.shape[0]

    guide = []; guideSet = set([])
    newD = []
    for d in range(nGuide):
        sel = -1
        minDist = 1e20
        dvec_bar = None
        dvec = D[d]

        # normalize
        scale = np.linalg.norm(dvec[-3:])
        dvec /= scale

        for hair in xrange(nHair):
            hvec = X[hair]
            diff = np.linalg.norm(dvec-hvec)
            if minDist > diff:
                minDist = diff
                sel = hair
                dvec_bar = dvec

        if sel not in guideSet:
            guide.append(sel)
            guideSet.add(sel)
            newD.append(dvec_bar)

    return list(guide), np.array(newD).transpose(), len(guide)


def guideSelect2016(fileName, nGuide):
    '''main function'''

    #debug
    nFrame = 5
    dump = DumpEngine("D:/tempDump")
    stage = 0
    setupDefaultLogger("D:/log/log.log")

    # global paramters
    lambda1 = para.lambda1
    xsima = para.xsima

    X, hairHeader, Data = SCGetMatrixAndHeader(fileName, nFrame) # X: len(u_s) x nHair

    if stage < 1:
        Us = np.asfortranarray(X, 'd')

        params = {'lambda1':lambda1, 'lambda2':0, 'return_model':True, 'model':None, 'posAlpha':True}
        D, ABi = spams.trainDL(Us, K=nGuide, iter=1, batchsize=10,**params) # D: len(u_s) x nGuide

        guide, D_bar, nGuide = pickGuideHair(D, X)
        params = {'lambda1':lambda1, 'lambda2':0, 'return_reg_path':False, 'pos':True}
        alpha = spams.lasso(Us, D = D_bar, **params) # alpha: nGuide x nHair
        alpha = alpha.transpose()

        # release half of the memory
        Us = None
        X = X.astype('f')

        guideSet = []
        for coef in alpha:
            tmp = []
            nnz = coef.nonzero()
            nnzCount = coef.count_nonzero()
            for i in range(nnzCount):
                val = coef.data[i]
                if val > xsima:
                    tmp.append(guide[nnz[1][i]])
            guideSet.append(tmp)


        # load all guide hair information into memory
        offset = (hairHeader.factor * 6, hairHeader.factor * 3, 3)
        nFrame = X.shape[0]/offset[0]
        # assert(nFrame == para.nFrame and hairHeader.factor == para.factor)

        dump.dump(1, (guide, offset, nFrame))
        dump.dump(2, guideSet)

    if stage < 2:
        if stage == 1:
            guide, offset, nFrame = dump.load(1)

        BgDictR = dict.fromkeys(guide) # nframe * (factor - 1) * 3 * 3
        BgDictT = dict.fromkeys(guide) # nframe * (factor - 1) * 3
        for g in guide:
            s0 = X[:offset[0],g].A1
            slotR = []
            slotT = []
            for i in range(nFrame):
                tmpR = []
                tmpT = []
                sp = X[offset[0]*i:offset[0]*(i+1), g].A1
                for j in range(1, hairHeader.factor):
                    ps0 = (s0[j*offset[2]:(j+1)*offset[2]], s0[offset[1]+j*offset[2]:offset[1]+(j+1)*offset[2]])
                    psp = (sp[j*offset[2]:(j+1)*offset[2]], sp[offset[1]+j*offset[2]:offset[1]+(j+1)*offset[2]])

                    tmpR.append(cd.vector_rotation_3D(ps0[1], psp[1]))
                    tmpT.append(psp[0] - ps0[0])
                slotR.append(tuple(tmpR))
                slotT.append(tuple(tmpT))
            BgDictR[g] = np.array(slotR)
            BgDictT[g] = np.array(slotT)

        dump.dump(3, (BgDictR, BgDictT))

    if stage < 3:
        if stage == 2:
            BgDictR, BgDictT= dump.load(3)
            guideSet = dump.load(2)
            guide, offset, nFrame = dump.load(1)

        # compute the weight
        cons = ({'type': 'eq',
                 'fun': lambda x: np.sum(x) - 1.0,
                 'jac': lambda x: np.ones(len(x))
                 },
                {'type': 'ineq',
                 'fun': lambda x: x,
                 'jac': lambda x: np.identity(len(x))
                 })

        sumError = 0.0
        for i in xrange(hairHeader.nHair):
            ghairs = guideSet[i]
            sfs = Data.getState(i).reshape((nFrame, 2, hairHeader.factor, 3))
            BgRs = np.array(map(lambda g: BgDictR[g], ghairs))
            BgTs = np.array(map(lambda g: BgDictT[g], ghairs))

            for j in range(1, hairHeader.factor): # we do not compute the weight of root particle
                sf = sfs[:, :, j, :]
                R = BgRs[:, :, j-1, :, :]
                T = BgTs[:, :, j-1, :]
                res, error = estimateWeight(sf, R, T, cons)
                print res.x, error
                sumError += error

if __name__ == "__main__":
    guideSelect2016(r"D:\Data\c0524\c0514.anim2", 200)