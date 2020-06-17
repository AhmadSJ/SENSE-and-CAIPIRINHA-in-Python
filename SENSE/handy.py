import numpy as np
import h5py
import matplotlib.pyplot as plt
import math

def caipi(r,ry,shifts):
    """Generates Caipirinha binary unit cells on sheared grids that form periodic lattices
    
    Parameters:
        r (int): total acceleration factor
        ry (int): acceleration factor in the y-direction
        shifts (list), (numpy array): size of shifts in the z-direction
    
    Returns:
        3D numpy array [shifts, ry, rz] containing the patterns
        
    
    """
    rz = int(r/ry)
    patterns = np.zeros((len(shifts),r,r))
    l = 0
    for i, shift in enumerate(shifts):
        cell = np.zeros((r,r))
        for row in range(rz):
            if (shift*row+rz)>= r:
                l = (shift*row%r+rz)%rz
            else:
                l = shift*row
            cell[row*ry,l::rz]=1
        patterns[i] = np.copy(cell)
    return patterns

def read_matlab(filename):
    """Load .mat file
    """
    def conv(path=''):
        p = path or '/'
        paths[p] = ret = {}
        for k, v in f[p].items():
            if type(v).__name__ == 'Group':
                ret[k] = conv(f'{path}/{k}')  # Nested struct
                continue
            v = v[()]  # It's a Numpy array now
            if v.dtype == 'object':
                # HDF5ObjectReferences are converted into a list of actual pointers
                ret[k] = [r and paths.get(f[r].name, f[r].name) for r in v.flat]
            else:
                # Matrices and other numeric arrays
                ret[k] = v if v.ndim < 2 else v.swapaxes(-1, -2)
        return ret

    paths = {}
    with h5py.File(filename, 'r') as f:
        return conv()


def plot_patterns(matrix,rows=8, cmap='gray'):
    n_patterns = len(matrix)
    R = matrix.shape[1]
    k = int(n_patterns/rows)+((n_patterns%rows)>0)
    s = np.arange(n_patterns)+1
    plt.figure(figsize=(12,10))
    for i in np.arange(n_patterns):
        plt.subplot(k,rows,s[i])
        plt.imshow(matrix[i],cmap=cmap)
        plt.title('Coil # {}'.format(i))
        ax = plt.gca();

        # Major ticks
        ax.set_xticks([]);
        ax.set_yticks([]);

        # Labels for major ticks
        ax.set_xticklabels([]);
        ax.set_yticklabels([]);

        # Minor ticks
        ax.set_xticks([], minor=True);
        ax.set_yticks([], minor=True);

        # Gridlines based on minor ticks
        ax.grid(which='minor', color='w', linestyle='-', linewidth=1)


def plot_caipi(matrix, shifts,rows=8):
    n_patterns = len(matrix)
    R = matrix.shape[1]
    k = int(n_patterns/rows)+((n_patterns%rows)>0)
    s = np.array(shifts)+1
    plt.figure(figsize=(16,8))
    for i in np.arange(n_patterns):
        plt.subplot(k,rows,s[i])
        plt.imshow(matrix[i])
        plt.title('$\Delta=$ {}'.format(shifts[i]))
        ax = plt.gca();

        # Major ticks
        ax.set_xticks(np.arange(0, R, 1));
        ax.set_yticks(np.arange(0, R, 1));

        # Labels for major ticks
        ax.set_xticklabels(np.arange(1, R+1, 1));
        ax.set_yticklabels(np.arange(1, R+1, 1));

        # Minor ticks
        ax.set_xticks(np.arange(-.5, R, 1), minor=True);
        ax.set_yticks(np.arange(-.5, R, 1), minor=True);

        # Gridlines based on minor ticks
        ax.grid(which='minor', color='w', linestyle='-', linewidth=1)

def swap(array,ind = 0):
    if array.ndim==3:
        if ind == 0:
            nc, x, y = array.shape
            h_array = np.zeros((x,y,nc),dtype=complex)
            for i in range(nc):
                h_array[:,:,i] = np.copy(array[i])
        else:
            x, y, nc = array.shape
            h_array = np.zeros((nc,x,y),dtype=complex)
            for i in range(nc):
                h_array[i] = np.copy(array[:,:,i])
    elif array.ndim==4:
        if ind == 0:
            nc, x, y, z = array.shape
            h_array = np.zeros((x,y,z,nc),dtype=complex)
            for i in range(nc):
                h_array[:,:,:,i] = np.copy(array[i])
        else:
            x, y, z, nc = array.shape
            h_array = np.zeros((nc,x,y,z),dtype=complex)
            for i in range(nc):
                h_array[i] = np.copy(array[:,:,:,i])
    return h_array

def take_center(array,size):
    """Takes the center data of an array

    Parameters:
        array (numpy array): array of data
        size (int): percentage of the array

    Returns (numpy array): center of array 
    """
    size /= 200
    Nc,Nx,Ny,Nz = array.shape
    nx,ny,nz = int(Nx/2), int(Ny/2),int(Nz/2)
    rx,ry,rz = math.ceil(Nx*size), math.ceil(Ny*size), math.ceil(Nz*size)
    reduced = array[:,nx-rx:nx+rx,ny-ry:ny+ry,nz-rz:nz+rz]
    #print(nx,ny,nz,rx,ry,rz)
    return reduced

def swappy(array):
    """Swap first and last dimensions of array
    """
    ndim = array.ndim
    count=1
    hulp = np.copy(array)
    for i in range(ndim-1)[::-1]:
        count+=1
        hulp = np.swapaxes(hulp,ndim-i-2,ndim-i-1)
    return hulp
