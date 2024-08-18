#pragma once

// -----------------------------------------------------------------------
// Collatz step function calculator 
// -----------------------------------------------------------------------

#include <stdint.h>
#include "gmp.h"

uint64_t fastest_collatz(mpz_t v);
