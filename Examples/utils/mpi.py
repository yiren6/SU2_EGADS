"""Module. Implements and configures the interface to MPI.

"""


try:
    from mpi4py import MPI
except ImportError:
    MPI = None

    COMM = None
    SIZE = 1
    RANK = 0
    ROOT = 0

else:
    COMM = MPI.COMM_WORLD
    SIZE = MPI.COMM_WORLD.Get_size()
    RANK = MPI.COMM_WORLD.Get_rank()
    ROOT = 0


__all__ = [
    'MPIPartitions',
    'COMM',
    'SIZE',
    'RANK',
    'ROOT',
    'isroot',
    'ifroot',
    'isrank',
    'ifrank',
    'isparallel',
    'ifparallel',
    'barrier',
    'broadcast',
    'allsum',
    'allgather',
    'gather',
    'scatter'
]


def isrank(rank):
    """Function. Return True on the root MPI rank, False elsewhere.

    """
    if rank is None:
        rank = ROOT

    return RANK == rank


def ifrank(function):
    """Function. Decorate a function to operate only on the given rank(s).

    """
    @functools.wraps(function)
    def conditional(*args, rank=None, **kwargs):
        """Internal function. Check if a function is run on a requested rank.

        """
        if isrank(rank):
            return function(*args, **kwargs)
        else:
            return None

    return conditional


def isroot():
    """Function. Return True on the root MPI rank, False elsewhere.

    """
    return RANK == ROOT


def ifroot(function):
    """Function. Decorate a function to operate only on the root rank.

    """
    @functools.wraps(function)
    def conditional(*args, **kwargs):
        """Internal function. Check if a function is run on the root rank.

        """
        if isroot():
            return function(*args, **kwargs)
        else:
            return None

    return conditional


def isparallel():
    """Function. Return True if more than one MPI rank is available.

    """
    return SIZE > 1


def ifparallel(function):
    """Function. Decorate a function to operate only on parallel runs.

    """
    @functools.wraps(function)
    def conditional(*args, **kwargs):
        """Internal function. Check if a function is run on a parallel run.

        """
        if isparallel():
            return function(*args, **kwargs)
        else:
            return None

    return conditional


def barrier():
    """Function. MPI synchronization barrier.

    """
    COMM.Barrier()


def broadcast(value):
    """Function. Broadcast a value from root to all ranks.

    """
    return COMM.bcast(value, root=ROOT)


def allsum(value):
    """Function. Sum a value over all ranks.

    """
    return COMM.allreduce(value, op=MPI.SUM)


def allgather(value):
    """Function. Gather a value over all ranks.

    """
    return COMM.allgather(value)


def gather(values):
    """Function. Gather partitioned values to root rank.

    """
    if not isparallel():
        return values

    return COMM.gather(values, root=ROOT)


def scatter(values):
    """Function. Scatter partitioned values to root rank.

    """
    if not isparallel():
        return values

    return COMM.scatter(values, root=ROOT)


