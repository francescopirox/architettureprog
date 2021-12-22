import array
import numpy as np
import scipy as sc;
import scipy.optimize;

d = 8; #dimensioni
nx = 64; #numero di pesci
iters = 250; #iterazioni

def dsload(filename, dtype = 'f'):
        with open(filename,'rb') as datafile:
                size = array.array('I');
                size.fromfile(datafile, 2);
                data = array.array(dtype);
                data.fromfile(datafile, size[0]*size[1]);   
        data = np.array(data);
        if(size[0]>1):
                return data.reshape(size[1],size[0]);
        return data;

def dssave(ds, filename, dtype = 'f'):
        if(ds.ndim == 2):
                n, d = ds.shape;
        else:
                n, d = len(ds), 1;
        with open(filename,'wb') as datafile:
                array.array('I', [d, n]).tofile(datafile);
                array.array(dtype, ds.reshape(np.prod(ds.shape))).tofile(datafile);

c = dsload('data/coeff32_'+str(d)+'.ds2');

def f(x): return np.exp(np.sum(x**2,axis=x.ndim-1))+ np.sum(x**2-c*x,axis=x.ndim-1);

xok = sc.optimize.minimize(f,-1/2*np.ones(d)).x;
print(xok);

dssave(xok.reshape(1,d),'sol32_'+str(d)+'_'+str(nx)+'_'+str(iters)+'.ds2',dtype='f');
dssave(xok.reshape(1,d),'sol64_'+str(d)+'_'+str(nx)+'_'+str(iters)+'.ds2',dtype='d');
