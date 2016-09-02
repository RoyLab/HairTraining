import spams
import numpy as np
from DataReader import HairDataReader, HairHeader
import local_para_small as para
import coordinates as cd
import sys
from scipy.optimize import minimize

import gc
import crash_on_ipy
import ipdb
import cPickle as pkl

def normalizePerColumn(X):
     X /= np.tile(np.sqrt((X*X).sum(axis=0)),(X.shape[0],1))


def evalError(x, m, lambda1):
    res = (((x * m[0]).dot(x) - 2 * m[1].dot(x) + m[2])*lambda1 + \
           (x * m[3]).dot(x) - 2 * m[4].dot(x) + m[5])
    return res

def evalDerive(x, m, lambda1):
    res = ((2 * m[0].dot(x) - 2 * m[1]).A1)*lambda1 + \
           (2 * m[0].dot(x) - 2 * m[1]).A1
    return res

def estimateWeight(sf, BgsfR, Bgsft, constrain):
    # sf: ndarray, nframe * ({pi}, {ti}), exclude 0, so n=factor-1
    # BgsfR/t: tuple, nframe * n * {R/t}
    lambda1 = 1.0

    nFrame = len(BgsfR)
    nGuide = len(BgsfR[0])

    pmat_square = np.matrix(np.zeros((nGuide, nGuide)))
    tmat_square = np.matrix(np.zeros((nGuide, nGuide)))
    plinearmat = np.zeros(nGuide)
    tlinearmat = np.zeros(nGuide)
    pcmat_square = 0.0
    tcmat_square = 0.0

    for i in range(1, nFrame):
        pc, tc = sf[i]
        Rs = BgsfR[i]
        Ts = Bgsft[i]
        pos_, dir_ = cd.point_trans(sf[0], Rs, Ts) # pos_, dir_ are flattened 1 x 3nGuide

        posMat = np.matrix(pos_.reshape(nGuide, 3))
        dirMat = np.matrix(dir_.reshape(nGuide, 3))
        pmat_square += posMat * posMat.T
        tmat_square += dirMat * dirMat.T
        plinearmat += pc[0] * posMat.T
        tlinearmat += tc[1] * dirMat.T
        pcmat_square += pc[0].dot(pc[0])
        tcmat_square += tc[0].dot(tc[0])

    res = minimize(evalError, [1.0 / nGuide] * nGuide, args=((pmat_square, plinearmat, pcmat_square,
            tmat_square, tlinearmat, tcmat_square), lambda1), jac=evalDerive, options={'disp': False},
                   method='SLSQP', constraints=constrain)




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
    return np.hstack(mat).transpose(), header

def pickGuideHair(D, X):
    return range(D.shape[1]), D


def guideSelect2016(fileName, nGuide):
    '''main function'''

    #debug
    nFrame = 2

    # global paramters
    lambda1 = para.lambda1
    xsima = para.xsima

    X, hairHeader = SCGetMatrixAndHeader(fileName, nFrame) # X: len(u_s) x nHair
    Us = np.asfortranarray(X, 'd')

    params = {'lambda1':lambda1, 'lambda2':0, 'return_model':True, 'model':None, 'posAlpha':True}
    D, ABi = spams.trainDL(Us, K=nGuide, iter=1, batchsize=10,**params) # D: len(u_s) x nGuide

    guide, D_bar = pickGuideHair(D, X)
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

    BgDictR = dict.fromkeys(guide)
    BgDictT = dict.fromkeys(guide)
    for g in guide:
        s0 = X[:offset[0],g]
        slotR = []
        slotT = []
        for i in range(nFrame):
            tmpR = []
            tmpT = []
            sp = X[offset[0]*i:offset[0]*(i+1), g].flatten()
            for j in range(1, hairHeader.factor):
                ps0 = (s0[j*offset[2]:(j+1)*offset[2]], s0[offset[1]+j*offset[2]:offset[1]+(j+1)*offset[2]])
                psp = (sp[j*offset[2]:(j+1)*offset[2]], sp[offset[1]+j*offset[2]:offset[1]+(j+1)*offset[2]])

                tmpR.append(cd.vector_rotation_3D(ps0[1], psp[1]))
                tmpT.append(psp[0] - ps0[0])
            slotR.append(tuple(tmpR))
            slotT.append(tuple(tmpT))
        BgDictR[g] = np.array(slotR)
        BgDictT[g] = np.array(slotT)

    # compute the weight
    cons = ({'type': 'eq',
             'fun': lambda x: np.sum(x) - 1.0,
             'jac': lambda x: np.ones(len(x))
             },
            {'type': 'ineq',
             'fun': lambda x: x,
             'jac': lambda x: np.identity(len(x))
             })

    for i in xrange(hairHeader.nHair):
        ghairs = guideSet[i]
        for j in range(1, hairHeader.factor): # we do not compute the weight of root particle
            estimateWeight(???????cons)





    pass




if __name__ == "__main__":
    guideSelect2016(r"D:\Data\c0524\c0514.anim2", 200)