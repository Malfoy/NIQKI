/*
 * This file comes more or less from the RedOak source code.
 */

#ifndef __COMMON_H__
#define __COMMON_H__

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <cstring>

// Tags for (debug) messages
#define MSG_LOG_HEADER \
  __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << ":"


// Make debug easy with macros:
// - DBG(instructions)
// - DEBUG_MSG(msg)
// - OMP_CRITICAL_DEBUG_MSG(msg)
#ifdef DEBUG
#  define DBG(instructions) instructions; (void)0
#  define DEBUG_MSG(msg)					\
  cout << MSG_LOG_HEADER << msg << endl;	\
  (void)0

#  ifdef NDEBUG
#    undef NDEBUG
#  endif


#else
#  define DBG(instructions)           (void)0
#  define DEBUG_MSG(msg)              (void)0
#endif

#endif
// Local Variables:
// mode:c++
// End:

