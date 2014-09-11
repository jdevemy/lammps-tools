# -*- coding: iso-8859-1 -*-
'''Python library with utils for script stuff'''

import itertools
import logging
import sys
import time
import os.path
import shelve
from multiprocessing import Pool
import multiprocessing.pool
try:
  from mpi4py import MPI
except ImportError:
  USE_MPI = False
else:
  USE_MPI = True

NB_MAX_CONFS = 10000000

RECURSION_LIMIT = 10000

pool = None

def init_multi(nb_proc):
  '''Init multi stuff'''

  # Up the recursion limit to pickle big mol (with atoms) in multiprocessing
  sys.setrecursionlimit(RECURSION_LIMIT)
  global pool
  pool = Pool(processes=nb_proc)
  logging.info('Launching %d processes', nb_proc)

def task_done(_):
  '''Callback when a parallel task has finished'''

  nb_tasks_waiting = pool._taskqueue.qsize()

  logging.debug('Task done : currently %d tasks waiting', nb_tasks_waiting)

def task_launch():
  '''A new parallel task has arrived'''

  nb_max_tasks_waiting = 2 * pool._processes

  nb_tasks_waiting = pool._taskqueue.qsize()

  logging.debug('Task launch : currently %d tasks waiting', nb_tasks_waiting)

  if nb_tasks_waiting > nb_max_tasks_waiting:
    logging.info('Too much tasks waiting (%d), waiting...', nb_tasks_waiting)
    while nb_tasks_waiting > nb_max_tasks_waiting:
      time.sleep(1)
      nb_tasks_waiting = pool._taskqueue.qsize()
      logging.debug('Waiting ... (%d)', nb_tasks_waiting)
    logging.info('Less tasks waiting (%d), restarting...', nb_tasks_waiting)

def close_multi():
  '''Close multi stuff'''

  if pool is not None:
    pool.terminate()

def init_mpi(func):
  '''Init MPI stuff, mainly for slaves'''

  if USE_MPI:
    comm = MPI.COMM_WORLD
    if comm.Get_size() > 1 and comm.Get_rank() != 0:
      rank = comm.Get_rank()
      logging.info('Slave %d: waiting', rank)
      # Main slave loop
      while True:
        params = comm.recv(source=0, tag=42)
        if params == None:
          logging.info('Slave %d: exiting', rank)
          MPI.Finalize()
          sys.exit(0)
        logging.info('Slave %d: start compute', rank)
        result = func(*params)
        logging.info('Slave %d: end compute', rank)
        comm.send((comm.Get_rank(), result), dest=0, tag=43)

def init_logging(verbosity):
  '''Init logging stuff'''

  logger = logging.getLogger()
  handler = logging.StreamHandler()
  formatter = logging.Formatter('%(asctime)s %(levelname)s [%(process)d]: %(message)s', '%d/%m/%Y-%H:%M:%S')
  handler.setFormatter(formatter)
  logger.addHandler(handler)

  if verbosity == 0:
    logger.setLevel(logging.WARNING)
  elif verbosity == 1:
    logger.setLevel(logging.INFO)
  else:
    logger.setLevel(logging.DEBUG)

def manage_arg_list(list_):
  '''Manage -l/--list arg'''

  step = 1
  i_cnfs = []

  # Step
  if '/' in list_:
    step = int(list_.split('/')[1])
    list_ = list_.split('/')[0]

  # Elems
  fields = list_.split(',')
  for elem in fields:
    if ':' in elem:
      if elem.startswith(':'):
        i_elem1 = 1
      else:
        i_elem1 = int(elem.split(':')[0])
      if elem.endswith(':'):
        i_elem2 = NB_MAX_CONFS
      else:
        i_elem2 = int(elem.split(':')[1])
      if i_elem1 <= 0 or i_elem2 <= 0:
        raise ValueError('Index should be positive')
      i_cnfs = itertools.chain(i_cnfs, xrange(i_elem1, i_elem2 + 1))
    else:
      i_elem = int(elem)
      if i_elem <= 0:
        raise ValueError('Index should be positive')
      i_cnfs = itertools.chain(i_cnfs, [i_elem])

  # Apply step
  if step != 1:
    i_cnfs = iter(list(i_cnfs)[::step])

  # Get the first i but reput it in the iterator
  i_first = i_cnfs.next()
  i_cnfs = itertools.chain([i_first], i_cnfs)

  return (i_first, i_cnfs)

def key2str(key):
  '''Convert list of object to a string key'''

  key_list = []

  # If a scalar, put it in a list for analyse
  if not isinstance(key, (list, tuple)):
    key = [key]

  for elem in key:
    if isinstance(elem, (str, unicode)):
      key_list.append(elem)
    elif isinstance(elem, (int, long, float, bool)):
      key_list.append(str(elem))
    elif isinstance(elem, (dict)):
      dkeys = elem.keys()
      dkeys.sort()
      for dkey in dkeys:
        key_list.append(dkey + ':' + key2str(elem[dkey]))
    elif isinstance(elem, (list, tuple)):
      for el in elem:
        key_list.append(key2str(el))
    elif isinstance(elem, (file)):
      mtime = os.path.getmtime(elem.name)
      key_list.append(elem.name + ':' + str(int(mtime)))
    else:
      key_list.append(str(elem))

  return '__'.join(key_list)

def compute(nb_proc, i_cnfs, dump, func, params, key=None):
  '''Main compute loop for conf analysis'''

  # Key is for persistent storage
  if key is not None:
    key_str = key2str(key)
    is_persist = True
  else:
    is_persist = False

  if is_persist:
    persist_filename = '.' + os.path.basename(sys.argv[0])
    persist = shelve.open(persist_filename, writeback=True)
    if not persist.has_key(key_str):
      logging.info('Creating persistent info to file %s', persist_filename)
      persist[key_str] = {}
    else:
      logging.info('Getting persistent info from file %s', persist_filename)
    logging.debug('Persist content for key %s : %s', key_str, persist[key_str])
    pers_data = persist[key_str]

  if nb_proc > 1:
    is_multi = True
  else:
    is_multi = False

  # Check MPI
  is_mpi = False
  if USE_MPI:
    comm = MPI.COMM_WORLD
    if comm.Get_size() > 1:
      is_mpi = True
      is_multi = False
      logging.info('Using MPI')
      mpi_slaves = [ None ] * (comm.Get_size() - 1)

  if is_multi:
    init_multi(nb_proc)

  results = []

  persist_results = {}

  nb_cnf = 0

  i_cnf = i_cnfs.next()
  # The dump is already sync to the first conf
  i_dump = i_cnf - 1

  # Loop on confs
  for cnf in dump:

    i_dump += 1

    # Not an interesting conf
    if i_dump != i_cnf:
      continue

    nb_cnf += 1

    if nb_cnf % 100 == 0:
      logging.info('Getting conf %d from file %s (%d\'th conf read)', i_cnf, cnf['filename'], nb_cnf)

    mtime = int(os.path.getmtime(cnf['filename']))
    cnf_key = (cnf['i'], cnf['filename'], mtime)

    # Main compute
    if is_persist and cnf_key in pers_data:
      logging.info('Result for conf %s get from persistent', cnf_key)
      result = pers_data[cnf_key]
    else:
      if is_multi:
        result = pool.apply_async(func, [cnf] + params, callback=task_done)
        # Store the waiting result to put it when ready to persist
        if is_persist:
          persist_results[cnf_key] = result
        task_launch()
      elif is_mpi:
        for i in xrange(len(mpi_slaves)):
          # Get the first free slave and send it the cnf
          if mpi_slaves[i] == None:
            mpi_slaves[i] = cnf_key
            comm.send([cnf] + params, dest=i + 1, tag=42)
            logging.info('Task %s sent to slave %d', i_cnf, i + 1)
            # Put a fake result
            result = 'SLAVE_%d' % (i + 1)
            break

        nb_free_slaves = mpi_slaves.count(None)
        logging.info('%d free slaves', nb_free_slaves)
        # If no more free slaves, waiting for (any) result
        if nb_free_slaves == 0:
          (rank, slave_result) = comm.recv(source=-1, tag=43)
          cnf_key = mpi_slaves[rank - 1]
          logging.info('Task %s received from slave %d', cnf_key, rank)

          # Special case for receiving the result of the previous send, don't use fake
          if result == 'SLAVE_%d' % rank:
            result = slave_result
          else:
            # Replace the fake result by a good one
            for i in xrange(len(results)):
              if results[i] == 'SLAVE_%d' % rank:
                results[i] = slave_result

          if is_persist:
            logging.info('Result for conf %s store on persistent', cnf_key)
            pers_data[cnf_key] = slave_result
            persist[key_str] = pers_data
            persist.sync()

          mpi_slaves[rank - 1] = None
      else:
        result = func(cnf, *params)
        if is_persist:
          logging.info('Result for conf %s store on persistent', cnf_key)
          pers_data[cnf_key] = result
          persist[key_str] = pers_data
          persist.sync()

    results.append(result)

    # Check if some procesors have failed with an exception
    if is_multi:
      for result in results:
        if isinstance(result, multiprocessing.pool.ApplyResult) and result.ready() and not result.successful():
          logging.critical('An exception occured on one conf, exit !')
          raise KeyboardInterrupt

    # Check if some persist parallel results are already ready
    if is_persist and is_multi:
      keys2del = []
      for (cnf_key, result) in persist_results.items():
        # If a result is ready, put it in persist
        if result.ready():
          logging.info('Result for conf %s store on persistent', cnf_key)
          pers_data[cnf_key] = result.get()
          persist[key_str] = pers_data
          persist.sync()
          keys2del.append(cnf_key)
      for cnf_key in keys2del:
        del persist_results[cnf_key]

    # Get i for next conf
    try:
      i_cnf = i_cnfs.next()
    except StopIteration:
      break

  if is_mpi:
    nb_free_slaves = mpi_slaves.count(None)
    # Wait for last slaves
    while nb_free_slaves != len(mpi_slaves):
      (rank, slave_result) = comm.recv(source=-1, tag=43)
      cnf_key = mpi_slaves[rank - 1]
      logging.info('Task %s received at end from slave %d', cnf_key, rank)
      # Replace the fake result by a good one
      for i in xrange(len(results)):
        if results[i] == 'SLAVE_%d' % rank:
          results[i] = slave_result
      if is_persist:
        logging.info('Result for conf %s store on persistent', cnf_key)
        pers_data[cnf_key] = slave_result
        persist[key_str] = pers_data
        persist.sync()

      mpi_slaves[rank - 1] = None
      nb_free_slaves = mpi_slaves.count(None)

    # Send shutdown message
    for i in xrange(len(mpi_slaves)):
      comm.send(None, dest=i + 1, tag=42)

    MPI.Finalize()

  # Get the results
  if is_multi:
    pool.close()

    logging.info('Waiting for the results')

    # Wait for the end of all tasks
    pool.join()

    # Convert results
    final_results = []
    for result in results:
      # Waiting result or direct ones (persist)
      if isinstance(result, multiprocessing.pool.ApplyResult):
        final_results.append(result.get())
      else:
        final_results.append(result)
    results = final_results

  # Store last parallel results in persist
  if is_persist and is_multi:
    keys2del = []
    for (cnf_key, result) in persist_results.items():
      # If a result is ready, put it in persist
      if result.ready():
        logging.info('Result for conf %s store on persistent', cnf_key)
        pers_data[cnf_key] = result.get()
        persist[key_str] = pers_data
        persist.sync()
        keys2del.append(cnf_key)
    for cnf_key in keys2del:
      del persist_results[cnf_key]

  if is_persist:
    logging.debug('Persist content for key %s : %s', key_str, persist[key_str])
    persist.close()

  logging.info('%d conf analysed', nb_cnf)

  return (nb_cnf, results)
