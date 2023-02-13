/*
 *    This file is part of the PowerBorn Library and SIMONA
 *    Please refer to the file LICENSE in either the PowerBorn package (LGPL)
 *    for SIMONA (proprietary) for License and Copyright information.
 */

#include <iostream>
#include <opencl_critical.h>

sig_atomic_t OpenClCriticalLock::my_signal = 0;

const char *OpenClCriticalLockException::what() const throw() { return "Error! Could not redirect all signals."; }

void OpenClCriticalLock::setSignal(int sig) { my_signal = sig; }

void OpenClCriticalLock::reset() {
  // std::cout << "destroying opencl critical lock" << std::endl;
  // restore previous state
  signal(SIGINT, sf[0]);
  signal(SIGSEGV, sf[1]);
  signal(SIGTERM, sf[2]);
  signal(SIGABRT, sf[3]);
#if _WIN32
  signal(SIGBREAK, sf[4]);
#else
  signal(SIGCHLD, sf[5]);
#endif
  // check if there was a signal, if yes, raise it
  // will propagate through nested locks outwards
}

void OpenClCriticalLock::release_check() {
  this->reset();
  if (my_signal == SIGSEGV || my_signal == SIGTERM || my_signal == SIGABRT) {
    std::cerr << "opencl critical section signal: " << my_signal << std::endl;
    throw OpenClCriticalLockSignal();
  }
}

OpenClCriticalLock::OpenClCriticalLock() {
  // std::cout << "creating opencl critical lock" << std::endl;
  for (unsigned int i = 0; i < 6; ++i) {
    sf[i] = SIG_DFL;
  }
  sf[0] = signal(SIGINT, OpenClCriticalLock::setSignal);
  sf[1] = signal(SIGSEGV, OpenClCriticalLock::setSignal);
  sf[2] = signal(SIGTERM, OpenClCriticalLock::setSignal);
  sf[3] = signal(SIGABRT, OpenClCriticalLock::setSignal);
#ifdef _WIN32
  sf[4] = signal(SIGBREAK, OpenClCriticalLock::setSignal);
#else
  sf[5] = signal(SIGCHLD, OpenClCriticalLock::setSignal);
#endif
  // test if signals will be successfully redirected
  for (unsigned int i = 0; i < 6; ++i) {
    if (sf[i] == SIG_ERR) {
      throw OpenClCriticalLockException();
    }
  }
}

OpenClCriticalLock::~OpenClCriticalLock() throw() { this->reset(); }
