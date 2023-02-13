/*
 *    This file is part of the PowerBorn Library and SIMONA
 *    Please refer to the file LICENSE in either the PowerBorn package (LGPL)
 *    for SIMONA (proprietary) for License and Copyright information.
 */

#ifndef OPENCL_CRITICAL_H_
#define OPENCL_CRITICAL_H_

#include <csignal>
#include <exception>

class OpenClCriticalLockException : public std::exception {
  const char *what() const throw();
};

class OpenClCriticalLockSignal : public std::exception {};

// this class is not thread safe, use only in one single thread!!!
// add static mutex and locks when ever my_signal is read or written to
class OpenClCriticalLock {
  typedef void (*signal_func)(int);

  static sig_atomic_t my_signal;
  signal_func sf[6];

  static void setSignal(int sig);
  void reset();

public:
  void release_check();
  OpenClCriticalLock();
  ~OpenClCriticalLock() throw();
};

#endif /* OPENCL_CRITICAL_H_ */
