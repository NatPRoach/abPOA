##  A header file to get you set going with Intel SIMD instrinsic programming.
##  <immintrin.h> is inlucded for SSE2, SSE41, AVX2 and AVX512F, AVX512BW
##  SSE4.1: floor and blend is available)
##  AVX2: double speed
##  AVX512F: quardruple speed
##  AVX512BW: byte and word operation

# Copyright (c) 2018 Mamy André-Ratsimbazafy
# Distributed under the Apache v2 License
# (license terms are at http://www.apache.org/licenses/LICENSE-2.0).
# This file may not be copied, modified, or distributed except according to
# those terms.

when defined(i386) or defined(amd64):
  # SIMD throughput and latency:
  #   - https://software.intel.com/sites/landingpage/IntrinsicsGuide/
  #   - https://www.agner.org/optimize/instruction_tables.pdf

  # Reminder: x86 is little-endian, order is [low part, high part]
  # Documentation at:
  #   https://software.intel.com/sites/landingpage/IntrinsicsGuide/

  when defined(vcc):
    {.pragma: x86_type, bycopy, header:"<intrin.h>".}
    {.pragma: x86, noDecl, header:"<intrin.h>".}
  else:
    {.pragma: x86_type, bycopy, header:"<x86intrin.h>".}
    {.pragma: x86, noDecl, header:"<x86intrin.h>".}
  type
    m128* {.importc: "__m128", x86_type.} = object
      raw: array[4, float32]
    m128d* {.importc: "__m128d", x86_type.} = object
      raw: array[2, float64]
    m128i* {.importc: "__m128i", x86_type.} = object
      raw: array[16, byte]
    m256* {.importc: "__m256", x86_type.} = object
      raw: array[8, float32]
    m256d* {.importc: "__m256d", x86_type.} = object
      raw: array[4, float64]
    m256i* {.importc: "__m256i", x86_type.} = object
      raw: array[32, byte]
    m512* {.importc: "__m512", x86_type.} = object
      raw: array[16, float32]
    m512d* {.importc: "__m512d", x86_type.} = object
      raw: array[8, float64]
    m512i* {.importc: "__m512i", x86_type.} = object
      raw: array[64, byte]
    mmask16* {.importc: "__mmask16", x86_type.} = distinct uint16
    mmask64* {.importc: "__mmask64", x86_type.} = distinct uint64
    # abPOA specific below here
    SIMDf* {.importc: "SIMDf", bycopy, header:"src/simd_instruction.h".} = m512
    SIMDi* {.importc: "SIMDi", bycopy, header:"src/simd_instruction.h".} = m512i
