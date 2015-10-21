
# coding: utf-8

# In[263]:

from image import MultiImage
if __name__ == "__main__":
    img = MultiImage("atlas/Christian/20150414-1-200-Tub-GFP/GFP-Controle_3/GFP-Controle_3_MMStack_Pos0_1.ome.tif")
    for i,frame in enumerate(img): 
        if i == 1:
            break
    #print frame.shapefrom image import MultiImage


# In[71]:

if __name__ == "__main__":

    a=hist(frame.flatten(),bins=200,cumulative=True,normed=True)

    for n,el in enumerate(a[0]):
        if el > 0.98:
            thres = a[1][n]
            print thres
            break
    print thres / (frame.max()-frame.min())
#a[0] > 0.95


# In[274]:

from numpy import *
from scipy import optimize

def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)
def lgaussianS(height, center_x, center_y):
    """Returns a gaussian function with the given parameters"""
    width_x = 2.
    width_y = 2.
    shift =0
    if height <= 0:
        shift=-height+0.01
    return lambda x,y: np.log(height+shift)                 -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2
def gaussianS(height, center_x, center_y):
    """Returns a gaussian function with the given parameters"""
    width_x = 2.
    width_y = 2.
    return lambda x,y: height*exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)
def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data, WDeux = True,log=False,cm=False):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    
    if cm:
        return np.array([params[0],params[1],params[2]]+[2,2])
   
    if WDeux:
        if log:
            datap = np.log(data)
            errorfunction = lambda p: ravel(lgaussianS(*p)(*indices(data.shape)) -
                                         datap )
            p = optimize.leastsq(errorfunction, params[:3],full_output =True)
        else:
            """
            errorfunction = lambda p: ravel(gaussianS(*p)(*indices(data.shape)) -
                                         data )
            p = optimize.leastsq(errorfunction, params[:3],full_output =True)"""
            #m = data.max()
            #data[data<0.30*m] = 0


            #indicesp = ma.array(data).nonzero()
            #xp,yp = indicesp
            xp,yp = indices(data.shape)
            errorfunction = lambda p: ravel(gaussianS(*p)(xp,yp) -
                                         data[xp,yp] )
            p = optimize.leastsq(errorfunction, params[:3],full_output =True)
        #print p[2]["nfev"], np.std(gaussianS(*p[0])(*indices(data.shape))-data)
        return np.array(p[0].tolist()+[2,2])
    else:
        errorfunction = lambda p: ravel(gaussian(*p)(*indices(data.shape)) -
                                     data)
        p, success = optimize.leastsq(errorfunction, params)
        return p


# In[271]:

from skimage.feature import blob_log,blob_dog,blob_doh
from skimage.feature import peak_local_max
from scipy import ndimage
from fast_blob_dog import blob_dog_fast
import time
def cm(frame):
    
    framep = np.array(frame,dtype=np.float)
    #print framep,
    T = np.sum(framep)
    cx = np.sum(np.arange(framep.shape[0],dtype=np.float)*np.sum(framep,axis=1))/T
    cy = np.sum(np.arange(framep.shape[1],dtype=np.float)*np.sum(framep,axis=0))/T
    #print cx,cy
    return cx,cy
def analyse_frame(frame,max_sigma=2,threshold=0.02,npix=4,timeit=False,cm=False,prune=False):
    if timeit : t0 = time.time()
    blobs =  blob_dog_fast(frame.copy(),max_sigma=max_sigma,threshold=threshold,overlap=0.01,prune=prune)#,num_sigma=1)
    if timeit : 
        print "Blob",time.time()-t0
        t0 = time.time()
    #ax.set_title()
    Xt,Yt = [],[]
    for blob in blobs:
        y0, x0, r = blob
        n=4
        if  n<x0<frame.shape[0]-n-1 and  n<y0<frame.shape[1]-n-1:
            #yp,xp =cm(frame[y0-n:y0+n+1,x0-n:x0+n+1])
            #print "log",
            params = fitgaussian(frame[y0-n:y0+n+1,x0-n:x0+n+1]-frame[y0-n:y0+n+1,x0-n:x0+n+1].min()+0.01,
                                 log=False,cm=cm)
            yp,xp=params[1:3]
            #n=5
            #print "Not log",
            #params2 = fitgaussian(frame[y0-n:y0+n+1,x0-n:x0+n+1],log=False)
            #yp1,xp1=params2[1:3]
            #print yp-yp1,xp-xp1
            #print yp,xp
            if not (0<xp<2*n and 0<yp<2*n):
                continue
            Xt.append(x0+xp-n)
            Yt.append(y0+yp-n)
            #print Xt[-1],Yt[-1]
            #figure()
            #imshow( frame[y0-n:y0+n+1,x0-n:x0+n+1], interpolation='none')
            #colorbar()
            #plot(xp,yp,"d",markersize=4)
            #print xp,yp
    #print "Out",time.time()-t0
    if timeit : 
        print "FitG",time.time()-t0
        t0 = time.time()
    return Xt,Yt
if __name__ == "__main__":

    f = figure(figsize=(10,10))
    #ax = f.add_subplot(121) 

    #ax.imshow(frame, interpolation='nearest')
    X,Y = analyse_frame(np.array(frame,dtype=np.float),max_sigma=2,threshold=1500,timeit=True)#,threshold=1000000)   

    print frame.dtype
    ax = f.add_subplot(111) 

    print len(X)
    ax.imshow(frame, interpolation='none')
    for x,y in zip(X,Y):
        # print x,y
        c = plt.Circle((x, y), 4, color="y", linewidth=2, fill=False)
        
        ax.add_patch(c)
        #p = plot(x,y,"o",markersize=2)
        #ax.add_patch(c)

    #ax.set_xlim(200,300)
    #ax.set_ylim(0,100)

    """
    for blob in peak_local_max(frame,threshold_abs=10000):
        y, x = blob
        #print x,y
        c = plt.Circle((x, y), 2, color="y", linewidth=2, fill=False)
    ax.add_patch(c)
"""
#imshow(frame)


# In[272]:

import multiprocessing as mp
from image import MultiImage
import PIL.Image as Image
import numpy as np
import types


def up_load_range(name,rrange,every=1,max_sigma=2,
                  threshold=0.02,multi=1,Zde=False,delta=0,lastc=-1,cm=False):
    
    #
    if delta == -1:
        #check the fit => fixed image
        if lastc != 0:
            rrange  = [ lastc+rrange[0]]
        #print "Check itnter , ", rrange
    if delta == 0:
        pass
    
    #print rrange,lastc
    if  type(name) == types.StringType:
        
        MI = MultiImage(name,conserve_memory=False)
        if len(MI)< len(rrange):
            rrange = range(len(MI))
        last = len(MI) 
    elif type(name) == types.ListType: 
        MI = [[] for x in name]
        
        if len(MI)< len(rrange):
            rrange = range(len(MI))
        timeit = False
        if timeit : 
            
            t0 = time.time()
        if multi == 1:
            #Else we load in multiproc
            for r in rrange:
                #print name[r]
                MI[r] = np.array(Image.open(name[r]),dtype=np.float)
        if timeit : 
            print "Load",time.time() -t0
            t0 = time.time()
            #print MI[r].shape,MI[r].dtype
        last = len(name)
        
    if delta > 0 and last-lastc < delta:
        #No update
        #print "No update"
        return [],[],[],[],last
    
    if delta <= 0 and lastc != 0:
        last=lastc
    
    #print rrange,lastc,delta,last
    if multi >1:
        
        output = mp.Queue()
        
        #lrange = range(-10-1,-1,1)
        #rrange=[lrange[:5],lrange[5:]]
        drange = np.array_split(np.array(rrange),multi)
        #print drange
        
        
        if type(name) == types.StringType:
            img = []
            for a in drange:
                init = rrange[0]
                img.append(np.zeros((len(a),MI[init].shape[0],MI[init].shape[1]),dtype=np.float))
                for count,i in enumerate(a):
                    img[-1][count] = MI[i].copy()
            processes = [mp.Process(target=load_range,
                            args=(img[x],
                                  max_sigma,threshold,output,cm,True)) for x in range(multi)]
            
        
        else:
            processes = [mp.Process(target=load_range,
                                args=([name[r] for r in drange[x]],
                                      max_sigma,threshold,output,cm,False)) for x in range(multi)]
        #processes = [mp.Process(target=testo,args=(name,rrange,max_sigma,threshold,output)) for x in range(1)]
        #processes = [mp.Process(target=testo, args=(5, output)) for x in range(4)]
        for p in processes:
            p.start()

            # Exit the completed processes
        #for p in processes:
        #    p.join()
            # Get process results from the output queue
       
        results = [output.get() for p in processes]
        X,Y,img = [],[],[]
        for Xt,Yt,imgp in results:
            X.extend(Xt)
            Y.extend(Yt)
            if imgp != []:
                img.append(imgp)
        #print np.array(np.concatenate(results,axis=0)).shape 
        #print results[0].shape
        #print results[0][0]
        return X,Y,[],np.mean(np.array(img),axis=0),last
  
    else:
        X,Y=[],[]
        Z=[]
        init = rrange[0]

        sumi = np.zeros_like(MI[init])
        for i in rrange:
            frame = MI[i]
            xtp,ytp = analyse_frame(np.array(frame,dtype=np.float),max_sigma=max_sigma,
                                    threshold=threshold,cm=cm,prune=False)
            X.extend(xtp)
            Y.extend(ytp)
            if Zde:
                Z.extend([i]*len(xtp))
            sumi+= frame
        #print len(X)
        if Zde:
            return X,Y ,Z,sumi/(1.0*len(rrange)),last
        return X,Y,[] ,sumi/(1.0*len(rrange)),last

def testo(**kwargs):
    print testo
def load_range(img,max_sigma=2,threshold=0.02,output=None,cm=False,MI=False):
    X,Y=[],[]
    imgf =[]
    #print img.shape
    
    for ni,frame in enumerate(img):
        
        #print frame.shape
        if not MI:
            frame = np.array(Image.open(frame),dtype=np.float)
        xtp,ytp = analyse_frame(frame,max_sigma=max_sigma,threshold=threshold,
                                timeit=False,cm=cm,prune=False)
        X.extend(xtp)
        Y.extend(ytp)  
        if ni == 0:
            imgf = frame.copy()
        else:
            imgf += frame
            
    #print "Put"
    if imgf != []:
        output.put([X,Y,imgf/len(img)])
    else:
        output.put([X,Y,imgf])
    
if __name__ == "__main__":
    import glob
    name = "../test-pil/atlas/Mickael/Francesca/09111_2A568_H3_Cy5_control_Cy5_1/"
    name = glob.glob(name+"/*.tif")
    name.sort()
    rangei = range(2960,2960+50,1) 
    #name = "../test-pil/atlas/Christian/20150414-1-200-Tub-GFP/GFP-Controle_3/GFP-Controle_3_MMStack_Pos0_1.ome.tif"
    #rangei=range(-5,-1,1)
    get_ipython().magic(u'timeit X,Y,Z,f,t = up_load_range(name,rrange=rangei,every=1,                                       max_sigma=2,threshold=500,multi=1,Zde=True,cm=False)')
    
    X,Y,Z,f,t = up_load_range(name,rrange=rangei,every=1,                                       max_sigma=2,threshold=500,multi=1,Zde=True,cm=False)
    print len(X)
    imshow(f)


# In[151]:

if __name__=="__main__":
    #print Z
    def ddist(delta,X,Y,Z):
        s = list(set(Z))
        s.sort()

        Z = np.array(Z)
        X= np.array(X)
        Y = np.array(Y)
        distance = []
        tos = s if delta == 0 else  s[:-delta]
        for el in tos:
            X1,Y1 = X[Z==el],Y[Z==el]
            X2,Y2 = X[Z==el+delta],Y[Z==el+delta]

            for x,y in zip(X1,Y1):
                distance.extend(np.array(np.sqrt((X2-x)**2+(Y2-y)**2)).tolist())
        distance = np.array(distance)

        if delta == 0:
            distance=distance[distance != 0]
        return distance*106
    D=[]
    for delta in [0,1,50,300]:#range(3):
        D.append(ddist(delta,X,Y,Z).copy())
        print len(D[-1])
        hist(D[-1],bins=100,range=[0,600],alpha=0.4)#,weights=np.ones_like(D[-1],dtype=np.float)/len(D[-1]))

    #for r in range


# In[153]:

class particle:
    def __init__(self,point,time):
        self.points=[]
        self.times =[]
        self.add(point,time)
    def __list__(self):
        return self.CM
    list = __list__
    def __len__(self):
        return len(self.times)
    def __getitem__(self, key):
        '''Vector elements can be accessed directly.'''
        return self.CM[key]
    def __repr__(self): 
        return str(len(self.times))+" "+str(self.CM)
    def add(self,point,time):
        self.points.append(point)
        self.times.append(time)
    
    @property
    def CM(self):
        return np.mean(np.array(self.points),axis=0)
if __name__ == "__main__":  
    a = np.array([particle([1,2],0)])
    a[0].add([0,0],1)
    print a,a[0][0]
    def filtering_particle(X,Y,Z):
        s = list(set(Z))
        s.sort()

        Z = np.array(Z)
        X= np.array(X)
        Y = np.array(Y)
        particle_list=[]
        old_p =[]
        tos = s
        el =tos[0]
        X1,Y1 = X[Z==el],Y[Z==el]
        for x,y in zip(X1,Y1):
            particle_list.append(particle([x,y],el))
        for el in tos:
            X2,Y2 = X[Z==el],Y[Z==el]
            for p in particle_list:
                if p.times[-1] in [el-1,el-2]:
                    #print "th"
                    D = 106*np.array(np.sqrt((X2-p[0])**2+(Y2-p[1])**2))
                    if np.sum(D < 100) ==1:

                        #print p, X2[D < 100] ,D[D<100]
                        p.add([X2[D < 100],Y2[D < 100]],el)

            #Remove lonely particle
            for p in particle_list:
                if p.times[-1] not in [el,el-1]:
                    if len(p) == 1:
                        particle_list.remove(p)
                    else:
                        old_p.append(p)
                        particle_list.remove(p)
        #filter the last frame
        for p in particle_list:
            if p.times[-1] != el:
                if len(p.times) == 1:
                    particle_list.remove(p)
        return old_p + particle_list
            

    P = filtering_particle(X,Y,Z)    
    print P[:10]
    A = np.array(P)
    N = [ len(p) for p in A]

    hist(N,bins=100)
    figure()
    print np.sum(N),len(X)
    plot([p[0] for p in A[::]],[p[1] for p in A[::]],"o")
        


# In[265]:

if __name__ == "__main__":
    from skimage.feature import peak_local_max
    from scipy import ndimage
    lbl = ndimage.label(frame)[0]
    
    a = np.array([[0,0,0,0,0],
                 [0,1,1,1,0],
                 [0,1,3,1,0],
                 [0,1,2,1,0],
                 [0,0,0,0,0]],dtype=np.float)
#    a = np.array([[ 11123. , 14880. , 13101.],
# [ 11944.,  18426. , 15685.],
# [ 11506.  ,14961.  ,14082.]],dtype=np.float)
    blobs =  blob_dog(a,max_sigma=2,threshold=0.01)
    def cm(framep):
        print framep
        T = np.sum(framep)
        cx = np.sum(np.arange(framep.shape[0],dtype=np.float)*np.sum(framep,axis=1))/T
        cy = np.sum(np.arange(framep.shape[1],dtype=np.float)*np.sum(framep,axis=0))/T
        return cx,cy
    print blobs
    y0,x0,r = blobs[0]
    #y0,x0=1,1
    
    n=2
    params = fitgaussian(a[y0-n:y0+n+1,x0-n:x0+n+1],log=True)
    yp,xp=params[1:3]
    x = x0+xp-2
    y = y0+yp-2
    print x,y
    imshow(a,interpolation='none')
    plot(x,y,"o",markersize=3)
    


# In[249]:

def test(N,nd,log=False,cm=False):
    dx = []
    dy = []
    X =[]
    Y =[]
    Ddelta = []
    fail = 0
    for i in range(N):
        cx,cy = np.random.rand(2)
        rescale = lambda x: 2*(1-2*x)
        cx = rescale(cx)
        cy= rescale(cy)
        test = gaussian(10000,nd+cx,nd+cy,2,2)(*indices((2*nd+1,2*nd+1))) +500*np.random.rand(2*nd+1,2*nd+1)
        test= test-test.min()+0.01
        #test = np.array(test,dtype=np.float)
        tos=False
        if tos:
            figure()
            imshow(test)
        #print test
        blobs =  blob_dog(test,max_sigma=2,threshold=100)
        y0,x0,r = blobs[0]
        n=4
        #print blobs
        #imshow(test[y0-n:y0+n+1,x0-n:x0+n+1])
        try:
            params = fitgaussian(test[y0-n:y0+n+1,x0-n:x0+n+1]-test[y0-n:y0+n+1,x0-n:x0+n+1].min()+0.01,
                                 log=log,cm=cm)
            delta = test[y0-n:y0+n+1,x0-n:x0+n+1]-test[y0-n:y0+n+1,x0-n:x0+n+1].min()+0.01 -                         gaussian(*params)(*indices((2*n+1,2*n+1))) 
            Ddelta.append(np.std(delta))
            yp,xp=params[1:3]

            x = x0+xp-n
            y = y0+yp-n

            #print "d",nd+cx-x,nd+cy-y
            dx.append(nd+cx-x)
            dy.append(nd+cy-y)
            X.append(x)
            Y.append(y)
            if tos:
                plot([x],[y],"o","g")
                plot([nd+cx,nd+cy],"o","r")

        except:
            fail += 1
    print fail
    return X,Y,dx ,dy,np.std(dy),np.mean(Ddelta)
    #%time 
    #%time test(1000,6,log=True)
if __name__ == "__main__":

    get_ipython().magic(u'time r=test(10000,10,cm=True)')
    get_ipython().magic(u'time r1=test(10000,10,log=False)')
    
    #plot(r[0])
    #hist(r[0])
    #hist(r[2],bins=20)
    


# In[252]:

if __name__ == "__main__":
    #print r
    a = hist(r1[0],bins=50)


# In[122]:

if __name__ == "__main__":
    """
    #a = "[[  6787.   8289.   7178.   7590.   6911.   6170.   6180.   7949.   5034.]
    # [  6809.   6229.   7810.   8214.   6346.   7751.   7047.   7933.   6234.]
    # [  7570.   7369.   7896.   9625.   9989.   7714.   8468.   6746.   6176.]
    # [  8100.   7775.   7719.  11123.  14880.  13101.   8201.   6816.   7552.]
    # [  8257.   7318.  10413.  11944.  18426.  15685.   9952.   6388.   7030.]
    # [  7751.   7124.   7422.  11506.  14961.  14082.   9550.   6859.   6825.]
    # [  6998.   8435.   7900.   7571.  12242.  10170.   8681.   6709.   6363.]
    ## [  8772.   6753.   6702.   8307.   7325.   7587.   6676.   6649.   6453.]
    # [  6755.   7365.   7514.   7367.   7866.   6537.   7208.   7205.   6749.]]"
    a = a.replace(".]","],")
    a=a.replace(".",",")
    a = a.replace(",]","]")
    #print a
    a = np.array([[  6787,   8289,   7178,   7590,   6911,   6170,   6180,   7949,   5034],
     [  6809,   6229,   7810,   8214,   6346,   7751,   7047,   7933,   6234],
     [  7570,   7369,   7896,   9625,   9989,   7714,   8468,   6746,   6176],
     [  8100,   7775,   7719,  11123,  14880,  13101,   8201,   6816,   7552],
     [  8257,   7318,  10413,  11944,  18426,  15685,   9952,   6388,   7030],
     [  7751,   7124,   7422,  11506,  14961,  14082,   9550,   6859,   6825],
     [  6998,   8435,   7900,   7571,  12242,  10170,   8681,   6709,   6363],
     [  8772,   6753,   6702,   8307,   7325,   7587,   6676,   6649,   6453],
     [  6755,   7365,   7514,   7367,   7866,   6537,   7208,   7205,   6749]],dtype=np.float)
    """
    params = fitgaussian(a)
    y,x=params[1:3]
    print params[1:3],cm(a)
    #x,y=cm(a)
    imshow(a,interpolation='none')
    plot(x,y,"o",markersize=3)
             


# In[115]:

if __name__ == "__main__":
    from sklearn.neighbors import KernelDensity
    xmin=0
    xmax=512
    ymin=0
    ymax=512
    Kd = KernelDensity(bandwidth=10.0,rtol=1e-4)
    X = np.array(X)
    Y = np.array(Y)
    print np.vstack([X,Y]).shape
    Kd.fit(np.array(np.vstack([X,Y]).T))
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    print positions.shape
    print Kd.score_samples(positions.T).shape
    f = np.reshape(Kd.score_samples(positions.T), xx.shape)


# In[116]:

if __name__ == "__main__":

    imshow(f.T,interpolation="none")
#scatter(X/5.12,Y/5.12)


# In[177]:

if __name__ == "__main__":

    res = 1
    npixel=512
    im = np.zeros((npixel*res,npixel*res))
    print len(X)
    fail = 0
    #X=np.array([5.5]*10)
    #Y=np.array([8]*10)
    X = np.array(X)
    Y = np.array(Y)

    for x,y in zip(X,Y):
        sigma=2.*res
        intensity=100.
        nd=5*res
        #print x,y
        try:
            GAU = gaussian(intensity,x-floor(x)+nd,y-floor(y)+nd,sigma,sigma)(*indices((2*nd+1,2*nd+1)))
            im[floor(x*res)-nd:floor(x*res)+nd+1,floor(y*res)-nd:floor(y*res)+nd+1] += GAU
        except:
            #print x,y
            fail += 1

    print fail


# In[174]:

if __name__ == "__main__":

    from matplotlib.colors import LogNorm

    figure(figsize=(20,20))

    imshow(im.T)#, norm=LogNorm())
    #scatter(X*res,Y*res)


# In[175]:

if __name__ == "__main__":

    plot(X,Y,"o")


# In[ ]:



