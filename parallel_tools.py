import os
import sys
import getpass
import logging
import traceback
import multiprocessing.pool

QUEUE_SIZE_MULTIPLIER = 8


class SyncAsyncPool:
  """A wrapper around multiprocessing.Pool which allows parallel processing but ordered results.
  This offers a compromise between synchronous and asynchronous processing, trying to get the
  benefits of both.
  It issues jobs asynchronously, but syncs up the results so they're processed in the same order
  the inputs were given. It does this by chunking the jobs, periodically stopping to wait for all
  jobs in the chunk to finish.
  It allows giving a callback which will be executed in the parent process. It will also be given
  results in the same order they were submitted."""

  def __init__(
    self, function, processes=None, queue_size=None, static_args=(), static_kwargs=None,
    callback=None, callback_args=()
  ):
    """Create a new SyncAsyncPool.
    processes can be None, "auto", an integer 0 or greater, or something that produces an integer
      >= 0 when converted with int(). None or "auto" mean to use as many subprocesses as there are
      cpu cores (reported by multiprocessing.cpu_count()). 0 means don't use any subprocesses. It
      will execute the function directly in the main process (and won't actually create a
      multiprocessing.Pool).
    queue_size can be None or an integer greater than 0. If it's None, the queue_size will be set
      to QUEUE_SIZE_MULTIPLIER * the number of processes."""
    # Validate arguments.
    if processes is not None and processes != 'auto':
      try:
        processes = int(processes)
      except (ValueError, TypeError):
        raise ValueError('processes must be an integer, None, or "auto" (received {!r})'.format(processes))
      if processes < 0:
        raise ValueError('processes must be 0 or greater (received {!r})'.format(processes))
    if processes == 'auto':
      processes = None
    if queue_size is not None and queue_size <= 0:
      raise ValueError('queue_size must be > 0 (received {!r})'.format(queue_size))
    # Are we actually doing multiprocessing, or should we do everything directly in one process?
    if processes == 0:
      self.multiproc = False
    else:
      self.multiproc = True
    if self.multiproc:
      self._pool = multiprocessing.pool.Pool(processes=processes)
    else:
      self._pool = None
    # Determine the number of processes.
    if processes is None or processes == 'auto':
      try:
        processes = multiprocessing.cpu_count()
      except NotImplementedError:
        processes = self._pool._processes
    self.processes = processes
    # Determine the queue size.
    if queue_size is None:
      if self.processes == 0:
        queue_size = 1
      else:
        queue_size = self.processes * QUEUE_SIZE_MULTIPLIER
    self.queue_size = queue_size
    self.function = function
    self.static_args = list(static_args)
    if static_kwargs is None:
      self.static_kwargs = {}
    else:
      self.static_kwargs = static_kwargs
    self.callback = callback
    self.callback_args = callback_args
    self.results = []

  def compute(self, *args, **kwargs):
    # Combine the static arguments with the args for this invocation.
    all_args = list(args) + self.static_args
    all_kwargs = self.static_kwargs.copy()
    all_kwargs.update(kwargs)
    # Send args to multiprocessing pool worker, or execute directly in this process if we're not
    # multiprocessing.
    if self.multiproc:
      result = self._pool.apply_async(with_context, [self.function]+all_args, all_kwargs)
    else:
      result = FakeResult(self.function(*all_args, **all_kwargs))
    self.results.append(result)
    if len(self.results) >= self.queue_size:
      self.flush()
      self.results = []

  def flush(self):
    if self.callback:
      for result in self.results:
        self.callback(result.get(), *self.callback_args)

  def close(self):
    if self.multiproc:
      self._pool.close()

  def join(self):
    if self.multiproc:
      self._pool.join()


class FakeResult(object):
  """A dummy version of multiprocessing.pool.AsyncResult to hold a result.
  It's convenient to have an object with the same API when we're not doing multiprocessing."""
  def __init__(self, result_data, timeout=None):
    self.result_data = result_data
  def get(self):
    return self.result_data


def with_context(fxn, *args, **kwargs):
  """Execute fxn, logging child process' stack trace for any Exceptions that are raised.
  When Exceptions are raised in a multiprocessing subprocess, the stack trace it gives ends where
  you call .get() on the .apply_async() return value.
  This logs the real stack trace and re-raises it.
  Usage:
  To execute real_fxn(arg1, arg2) through this, do:
    result = pool.apply_async(with_context, args=(real_fxn, arg1, arg2))
  Note: This would be a decorator, but that doesn't work with multiprocessing.
  Functions below the top-level of a module aren't picklable and so can't be passed through
  a Queue to subprocesses:
  https://stackoverflow.com/questions/8804830/python-multiprocessing-pickling-error/8805244#8805244
  NOTE: This must execute in the child process.
  """
  try:
    return fxn(*args, **kwargs)
  except Exception as exception:
    tb = traceback.format_exc()
    logging.critical(tb)
    exception.child_context = get_exception_data()
    raise exception


def get_exception_data():
  exception_type, exception_value, trace = sys.exc_info()
  trace_events = []
  for filename, lineno, fxn_name, code in traceback.extract_tb(trace):
    entry = {'file':filename, 'line':lineno, 'function':fxn_name, 'code':code}
    trace_events.append(entry)
  del trace
  return {'type':exception_type.__name__, 'traceback':trace_events}


def format_traceback(exception_data):
  lines = ['Traceback (most recent call last):']
  for entry in exception_data['traceback']:
    lines.append('  File "{file}", line {line}, in {function}'.format(**entry))
    lines.append('    '+entry['code'])
  return '\n'.join(lines)


def scrub_tb_paths(exception_data, script_path=__file__):
  """Convenience function that applies scrub_paths() to the paths in the traceback of data returned
  from get_exception_data()."""
  trace_events = exception_data['traceback']
  paths = scrub_paths([entry['file'] for entry in trace_events], script_path=script_path)
  for path, entry in zip(paths, trace_events):
    entry['file'] = path
  return exception_data


def scrub_paths(paths, script_path=__file__):
  """Do your best to remove sensitive data from paths."""
  script_dirs = get_script_dirs(script_path=script_path)
  # First, remove the common start of the paths.
  paths = abbreviate_paths(paths, keep_last=True)
  for path in paths:
    parts = path.split(os.sep)
    # If the start of the path matches our script_dir, just chop that off and return it.
    matched = False
    for script_dir in script_dirs.values():
      if path.startswith(script_dir):
        tail = path[len(script_dir)+1:]
        # But keep our immediate directory name.
        dirname = script_dir.split(os.sep)[-1]
        path = os.path.join(dirname, tail)
        matched = True
        break
    if matched:
      yield path
      continue
    # If the current user's username appears in the path, remove it and everything before it.
    username = getpass.getuser()
    for i, part in enumerate(parts):
      if part == username:
        parts = parts[i+1:]
        break
    # If the path starts with "/home/[something]" or "/user/[something]", remove those parts.
    if len(parts) >= 4 and (parts[:2] == ['', 'home'] or parts[:2] == ['', 'user']):
      parts = parts[3:]
    # Finally, keep only the last 2 directory names (unless the result would be longer than what
    # the previous steps produced).
    if len(parts) > 3:
      parts = path.split(os.sep)
      parts = parts[len(parts)-3:]
    yield os.sep.join(parts)


def get_script_dirs(script_path=__file__):
  """Try to get every possible representation of this script's directory.
  Try every variation of absolute, relative, following symbolic links and not."""
  script_dirs = {}
  script_dirs['rel_literal'] = os.path.dirname(script_path)
  script_dirs['abs_literal'] = os.path.abspath(script_dirs['rel_literal'])
  script_dirs['abs_linked'] = os.path.dirname(os.path.realpath(script_path))
  script_dirs['rel_linked'] = os.path.relpath(script_dirs['abs_linked'])
  return script_dirs


def abbreviate_paths(paths, keep_last=False):
  """Shorten paths by removing their common prefix.
  If keep_last, then never remove the last path component (which would happen if all the paths were
  the same)."""
  longest_prefix = get_longest_path_prefix(paths, return_type='list')
  for path in paths:
    parts = path.split(os.sep)
    if longest_prefix:
      starting_i = len(longest_prefix)
      if starting_i == len(parts) and keep_last:
        starting_i -= 1
      parts = parts[starting_i:]
      yield os.sep.join(parts)
    else:
      yield path


def get_longest_path_prefix(paths, return_type='list'):
  """Get the longest shared starting path among a set of paths.
  /home/me/.local/thing.txt, /home/me/.local/other.txt, /home/me/.cache/asdf.txt -> /home/me
  path/to/place/, path/to/other/place/with/thing.txt, path/to/here.txt -> path/to
  one/two/three, four/five/six.txt, four/seven/eight -> ''
  Mixing absolute and relative paths isn't an exception, but it always returns ''.
  If return_type == 'list', return the list of path components.
  If return_type == 'str', return the components joined on os.sep."""
  longest_prefix = None
  for path in paths:
    parts = path.split(os.sep)
    if longest_prefix is None:
      longest_prefix = parts
    else:
      new_longest = []
      for part1, part2 in zip(parts, longest_prefix):
        if part1 == part2:
          new_longest.append(part1)
      longest_prefix = new_longest
  if return_type == 'list':
    return longest_prefix
  else:
    return os.sep.join(longest_prefix)
