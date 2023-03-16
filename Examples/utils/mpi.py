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


class MPIPartitions:
    """Dataclass. MPI data partitioning information container (for all ranks).

    """
    def __init__(self, ids, domain):
        # distinguish between internal (domain) and exchange (halo) indices
        internal_local = []
        exchange_local = []

        for index, flag in enumerate(domain):
            if flag:
                internal_local.append(index)
            else:
                exchange_local.append(index)

        if isparallel():
            # gather sub-sets of IDs known to each rank to all ranks
            distribution = COMM.allgather(ids)

            # gather internal indices to all ranks
            internals_local = COMM.allgather(internal_local)
            internals_total = []

            # gather exchange indices to all ranks
            exchanges_local = COMM.allgather(exchange_local)
            exchanges_total = []

            # gather unique internal IDs found across all rank to all ranks
            gathered = []

            for ids, internal_local, exchange_local in zip(distribution, internals_local, exchanges_local):
                for index in internal_local:
                    gathered.append(ids[index])

            # calculate total internal and exchange indices for each rank
            for ids, internal_local, exchange_local in zip(distribution, internals_local, exchanges_local):

                # get internal and exchange IDs
                internal_ids = [
                    ids[local] for local in internal_local
                ]
                exchange_ids = [
                    ids[local] for local in exchange_local
                ]

                # match internal IDs to gathered IDs
                internal_total = []
                for local, _id in zip(internal_local, internal_ids):
                    total = gathered.index(_id)

                    internal_total.append(total)

                # match exchange IDs to gathered IDs
                exchange_total = []
                for local, _id in zip(exchange_local, exchange_ids):
                    total = gathered.index(_id)

                    exchange_total.append(total)

                internals_total.append(internal_total)
                exchanges_total.append(exchange_total)

        else:
            distribution = [
                ids
            ]
            gathered = ids

            internals_local = [
                internal_local
            ]
            internals_total = [
                internal_local
            ]
            exchanges_local = [
                exchange_local  # should be empty
            ]
            exchanges_total = [
                exchange_local  # should be empty
            ]

        # map local indices to total indices
        internals = [
            dict(zip(local, total)) for (local, total) in zip(internals_local, internals_total)
        ]
        exchanges = [
            dict(zip(local, total)) for (local, total) in zip(exchanges_local, exchanges_total)
        ]

        # store partitioned IDs (only on root)
        self._ids = list(gathered)

        # store partitioning data for each rank
        self._partitions = []

        for rank, (ids, internal, exchange) in enumerate(zip(distribution, internals, exchanges)):
            partition = MPIPartition(rank, ids, internal, exchange)

            self._partitions.append(partition)

        # total size of the partitioned data
        self._size = len(gathered)

    def __repr__(self):
        # allow string representation
        return f'{self.__class__.__name__}({self.size}, partitions={pprint.pformat(self._partitions)})'

    def __iter__(self):
        # allow iteration
        return iter(self._partitions)

    def __len__(self):
        # allow length operator
        return len(self._partitions)

    def __getitem__(self, rank):
        # allow index access
        try:
            return self._partitions[rank]
        except IndexError as err:
            raise IndexError(f'Rank: {rank} exceeds size {SIZE}') from err

    @property
    def size(self):
        """Property. Combined size of internal and exchange vertices.

        """
        return self._size

    @property
    def current(self):
        """Property. Partition of the current rank.

        """
        return self._partitions[RANK]

    @property
    def ids(self):
        """Property. List of IDs owned by this partition.

        """
        return self._ids

    def gather(self, array):
        """Function. Gather partitioned arrays to root rank.

        """
        array = np.ascontiguousarray(array)

        if not isparallel():
            return array  # short-circuit the function on serial runs

        dtype = array.dtype
        shape = array.shape

        # allocate send (local) buffer
        send = array[self.current.internal.local]

        # allocate recv (total) buffer
        if isroot():
            size = self.size
        else:
            size = 0

        _, *ndim = shape

        recv = np.empty((size, *ndim), dtype=dtype)

        # non-blocking send the data on all ranks ...
        request = COMM.Isend(send, dest=ROOT, tag=0)

        # ... blocking receive the data on root (from each rank)
        if isroot():
            for partition in self._partitions:
                temp = np.empty((partition.internal.size, *ndim), dtype=dtype)

                COMM.Recv(temp, source=partition.rank, tag=0)

                recv[partition.internal.total] = temp

        request.wait()

        return recv

    def scatter(self, array):
        """Function. Scatter values to rank partitions.

        """
        send = np.ascontiguousarray(array)

        if not isparallel():
            return send  # short-circuit the function on serial runs

        # calculate data properties and broadcast them to all ranks
        if isroot():
            dtype = send.dtype
            shape = send.shape
        else:
            dtype = None
            shape = None

        dtype = broadcast(dtype)
        shape = broadcast(shape)

        size, *ndim = shape

        if not size == self.size:
            raise ValueError(f"Expected size {self.size}, got {size}")

        # allocate recv (local) buffer
        recv = np.empty((self.current.size, *ndim), dtype=dtype)
        recv_internal = np.empty((self.current.internal.size, *ndim), dtype=dtype)
        recv_exchange = np.empty((self.current.exchange.size, *ndim), dtype=dtype)

        # non-blocking recv the data on all ranks ...
        request_internal = COMM.Irecv(recv_internal, source=ROOT, tag=1)
        request_exchange = COMM.Irecv(recv_exchange, source=ROOT, tag=2)

        # ... blocking send the data on root
        if isroot():
            for partition in self._partitions:
                send_internal = np.ascontiguousarray(send[partition.internal.total])
                send_exchange = np.ascontiguousarray(send[partition.exchange.total])

                COMM.Send(send_internal, dest=partition.rank, tag=1)
                COMM.Send(send_exchange, dest=partition.rank, tag=2)

        request_internal.wait()
        request_exchange.wait()

        # return the receiving buffer (now hopefully filled with data)
        recv[self.current.internal.local] = recv_internal
        recv[self.current.exchange.local] = recv_exchange

        return recv

    def reorder(self, ids):
        """Method. Re-order the partitions based on new ordering of IDs.

        """
        ordered_distribution = COMM.allgather(ids)
        current_distribution = [
            partition.ids for partition in self
        ]

        # loop once to compute new ordering of IDs in total array
        current_gathered = []
        ordered_gathered = []
        for partition, current_ids, ordered_ids in zip(self, current_distribution, ordered_distribution):
            for index in partition.internal.local:
                current_gathered.append(current_ids[index])
                ordered_gathered.append(ordered_ids[index])

        total = [
            current_gathered.index(_id) for _id in ordered_gathered
        ]

        # loop once to compute new ordering of IDs in local array
        for partition, current_ids, ordered_ids in zip(self, current_distribution, ordered_distribution):
            local = [
                current_ids.index(_id) for _id in ordered_ids
            ]

            partition.reorder(local, total)

        # re-order IDs
        self._ids = [
            self._ids[index] for index in total
        ]

        return total


class MPIPartition:
    """Dataclass. MPI data partitioning information (for one rank).

    """
    def __init__(self, rank, ids, internal, exchange=None):
        # rank for which partitioning data is being stored
        self._rank = rank

        # mapping of local internal (or "physical", or "domain") indices to total indices
        self._internal = MPIPartitionView(internal)
        # mapping of local exchange (or "ghost", or "halo") indices to total indices
        self._exchange = MPIPartitionView(exchange or {})

        # list of ids owned by this partition
        self._ids = list(ids)

    def __repr__(self):
        # allow string representation
        if self.exchange.size == 0:
            return f'{self.__class__.__name__}(internal={self.internal.size})'
        else:
            return f'{self.__class__.__name__}(internal={self.internal.size}, exchange={self.exchange.size})'

    @property
    def rank(self):
        """Property. Rank for which partitioning data was calculated.

        """
        return self._rank

    @property
    def size(self):
        """Property. Combined size of internal and exchange vertices.

        """
        return self.internal.size + self.exchange.size

    @property
    def internal(self):
        """Property. Mapping of internal (i.e. physical, domain) indices.

        """
        return self._internal

    @property
    def exchange(self):
        """Property. Mapping of exchange (i.e. ghost, halo) indices.

        """
        return self._exchange

    @property
    def ids(self):
        """Property. List of IDs owned by this partition.

        """
        return self._ids

    def reorder(self, local, total):
        """Method. Re-order the IDs and index maps.

        """
        self._ids = [
            self._ids[index] for index in local
        ]

        self.internal.reorder(local, total)
        self.exchange.reorder(local, total)


class MPIPartitionView:
    """Dataclass. MPI partition index map (local and global).

    """
    def __init__(self, mapping):
        # length of the partition subset
        size = len(mapping)

        self._size = size

        # local and total indices
        if size == 0:
            local, total = [], []
        else:
            local, total = zip(*mapping.items())

        self._local = np.asarray(local, dtype=int)
        self._total = np.asarray(total, dtype=int)

    @property
    def size(self):
        """Property. Return the size of the partition mapping.

        """
        return self._size

    @property
    def local(self):
        """Property. Return the local indices of the partition mapping.

        """
        return self._local

    @property
    def total(self):
        """Property. Return the total indices of the partition mapping.

        """
        return self._total

    def reorder(self, local, total):
        """Method. Re-order the local and total indices.

        """
        order = [
            total[index] for index in self._total
        ]

        self._total = np.array(order, dtype=int)

        order = [
            local[index] for index in self._local
        ]

        self._total = np.array(order, dtype=int)
