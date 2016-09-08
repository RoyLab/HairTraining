import spams
import numpy as np
from DataReader import HairDataReader, HairHeader
import local_para_small as para
import coordinates as cd
import sys
from scipy.optimize import minimize
from common_tools import *
from progressbar import ProgressBar
from multiprocessing import Pool, Manager, Lock

import crash_on_ipy

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

    #print nFrame, nGuide

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
                   options={'disp': False, 'ftol':1e-12, 'maxiter':100}, method='SLSQP', constraints=constrain, )

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


def runSCG(args):
    fileName, spcereg, offset, seq, batch, l = args
    reader = HairDataReader(fileName, {'type':'anim2'})
    mat = []
    reader.seek(seq*batch)
    for i in xrange(batch):
        frame = reader.getNextFrameNoRewind()
        if not frame: break

        pos = np.array(frame.position)
        dir = np.array(frame.direction)
        rigid = frame.headMotion
        invR = np.linalg.inv(rigid[0])
        pos, dir = cd.inverseRigidTrans(invR, rigid[1], pos, dir, batch=True)
        pos.shape = -1, offset
        dir.shape = -1, offset
        mat.append(pos*spcereg)
        mat.append(dir)

    reader.close()
    printMP(l, "Finish read batch %d..." % seq)
    return seq, mat


def SCGetMatrixAndHeaderMP(fileName, nFrame=None):
    reader = HairDataReader(fileName, {'type':'anim2'})
    factor = para.factor
    offset = factor * 3
    spcereg = para.lambda_balance  # regularize the space item
    if not nFrame:
        nFrame = reader.nFrame

    pool = Pool()
    m = Manager()
    l = m.Lock()
    batchsz = 10
    job_args = [(fileName, spcereg, offset, i, batchsz, l) for i in xrange(nFrame/batchsz) ]
    mat = pool.map(runSCG, job_args)

    pool.close()
    pool.join()
    reader.close()

    mat.sort(key=lambda x: x[0])
    mat2 = []
    for item in mat:
        for arr in item[1]:
            mat2.append(arr)

    mat = None
    X = np.hstack(mat2).transpose()
    mat2 = None

    print "\rFinished Reading!"

    header = HairHeader()
    header.nParticle = reader.nParticle
    header.factor = factor
    header.nHair = header.nParticle / factor

    return X, header, XWrapper(X, offset*2)

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

    #return range(nGuide), D.transpose(), nGuide

    guide = []; guideSet = set([])
    newD = []
    print "Guide hair selection..."
    bar = ProgressBar().start()
    for d in range(nGuide):
        sel = -1
        minDist = 1e20
        dvec_bar = None
        dvec = D[d]

        # normalize
        scale = np.linalg.norm(dvec[-3:])
        dvec /= scale

        for hair in xrange(nHair):
            hvec = X[hair].A1
            diff = np.linalg.norm(dvec-hvec)
            if minDist > diff:
                minDist = diff
                sel = hair
                dvec_bar = hvec

        if sel not in guideSet:
            guide.append(sel)
            guideSet.add(sel)
            newD.append(dvec_bar)
        bar.update(100*(d+1)/nGuide)

    bar.finish()

    return np.array(guide), np.asfortranarray(np.array(newD, 'd').transpose()), len(guide)

def genMatrixFromGuide(guide, X):
    newD = []
    X = X.transpose()
    for id in guide:
        dvec = X[id].A1
        newD.append(dvec)

    return np.asfortranarray(np.array(newD, 'd').transpose())

def selectByRandom(nGuide, X):
    from common_tools import genRandomNonRepeatArray
    N = X.shape[1]
    arr = genRandomNonRepeatArray(nGuide, N)

    return arr, genMatrixFromGuide(arr, X), nGuide

def selectByChai2016(nGuide, X, initD=None):
    lambda1 = para.lambda1
    Us = np.asfortranarray(X, 'd')

    params = {'lambda1': lambda1, 'lambda2': 0, 'return_model': True, 'model': None, 'posAlpha': True}
    D, ABi = spams.trainDL(Us, D=initD, K=nGuide, iter=100, batchsize=10, **params)  # D: len(u_s) x nGuide

    guide, D_bar, nGuide = pickGuideHair(D, X)

    print "Got %d guide hairs" % nGuide
    return guide, D_bar, nGuide

def guideSelect2016(fileName, nGuide):
    '''main function'''

    #debug
    nFrame = 200
    dump = DumpEngine("D:/tempDump")
    stage = 2
    setupDefaultLogger("D:/log/log.log")

    # global paramters
    xsima = para.xsima
    lambda1 = para.lambda1

    X, hairHeader, Data = SCGetMatrixAndHeader(fileName, nFrame) # X: len(u_s) x nHair
    X = X.astype('f')

    if stage < 1:
        startTime = getTimeStr()

        guide, D_bar, nGuide = selectByRandom(nGuide, X)

        params = {'lambda1':lambda1, 'lambda2':0, 'return_reg_path':False, 'pos':True}
        Us = np.asfortranarray(X, 'd')
        alpha = spams.lasso(Us, D = D_bar, **params) # alpha: nGuide x nHair
        alpha = alpha.transpose()

        # release half of the memory
        Us = None

        guideSet = []
        print "Guide set per normal hair..."
        bar = ProgressBar().start()
        count = 0
        for coef in alpha:
            count += 1
            tmp = []
            nnz = coef.nonzero()
            nnzCount = coef.count_nonzero()
            for i in range(nnzCount):
                val = coef.data[i]
                if val > xsima:
                    tmp.append(guide[nnz[1][i]])
            guideSet.append(tmp)
            bar.update(100*count/hairHeader.nHair)
        bar.finish()


        # load all guide hair information into memory
        offset = (hairHeader.factor * 6, hairHeader.factor * 3, 3)
        nFrame = X.shape[0]/offset[0]
        # assert(nFrame == para.nFrame and hairHeader.factor == para.factor)

        dump.dump(1, (guide, offset, nFrame))
        dump.dump(2, guideSet)

        print "Stage 0: start from "+ startTime
        print "Stage 0: ended at "+getTimeStr()

        exit(0)

    if stage < 2:
        if stage == 1:
            guide, offset, nFrame = dump.load(1)

        BgDictR = dict.fromkeys(guide) # nframe * (factor - 1) * 3 * 3
        BgDictT = dict.fromkeys(guide) # nframe * (factor - 1) * 3
        guideList = guide.tolist()
        for g in guide:
            print "load guide hair matrix %d ..." % guideList.index(g)
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
            if i in guide:
                continue

            print "Compute weight of normal %d ..." % i

            ghairs = guideSet[i]
            sfs = Data.getState(i).reshape((nFrame, 2, hairHeader.factor, 3))
            BgRs = np.array(map(lambda g: BgDictR[g], ghairs))
            BgTs = np.array(map(lambda g: BgDictT[g], ghairs))

            for j in range(1, hairHeader.factor): # we do not compute the weight of root particle
                sf = sfs[:, :, j, :]
                R = BgRs[:, :, j-1, :, :]
                T = BgTs[:, :, j-1, :]
                res, error = estimateWeight(sf, R, T, cons)
                sumError += error

        print "sum of error %f" % sumError

if __name__ == "__main__":
    nFrame = 200
    fileName = r"D:\Data\c0524\c0514.anim2"
    # X, hairHeader, Data = SCGetMatrixAndHeader(fileName, nFrame)  # X: len(u_s) x nHair

    import time
    t = time.time()
    SCGetMatrixAndHeaderMP(fileName, nFrame)  # X: len(u_s) x nHair
    print "MP: %f" % (time.time() - t)

    t = time.time()
    SCGetMatrixAndHeader(fileName, nFrame)
    print "SP: %f" % (time.time() - t)
