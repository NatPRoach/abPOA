import simd_instruction

import ../../conduit/src/poGraphUtils

const
  ABPOA_GLOBAL_MODE* = 0
  ABPOA_LOCAL_MODE* = 1
  ABPOA_EXTEND_MODE* = 2

## #define ABPOA_SEMI_MODE 3
##  gap mode

const
  ABPOA_LINEAR_GAP* = 0
  ABPOA_AFFINE_GAP* = 1
  ABPOA_CONVEX_GAP* = 2
  ABPOA_EXTRA_B* = 10
  ABPOA_EXTRA_F* = 0.01
  ABPOA_CIGAR_STR* = "MIDXSH"
  ABPOA_CMATCH* = 0
  ABPOA_CINS* = 1
  ABPOA_CDEL* = 2
  ABPOA_CDIFF* = 3
  ABPOA_CSOFT_CLIP* = 4
  ABPOA_CHARD_CLIP* = 5
  ABPOA_SRC_NODE_ID* = 0
  ABPOA_SINK_NODE_ID* = 1
  ABPOA_OUT_CONS* = 0
  ABPOA_OUT_MSA* = 1
  ABPOA_OUT_CONS_MSA* = 2
  ABPOA_OUT_GFA* = 3
  ABPOA_OUT_CONS_GFA* = 4
  ABPOA_HB* = 0
  ABPOA_HC* = 1
  ABPOA_MF* = 2

##  NOTE: upper boundary of in_edge_n is pow(2,30)
##  for MATCH/MISMATCH: node_id << 34  | query_id << 4 | op
##  for INSERTION:      query_id << 34 | op_len << 4   | op
##  for DELETION:       node_id << 34  | op_len << 4   | op 
##    op_len is always equal to 1
##  for CLIP            query_id << 34 | op_len << 4   | op

type
  uint8_t  {.header: "<stdint.h>" importc: "uint8_t".} = uint8
  uint16_t {.header: "<stdint.h>" importc: "uint16_t".} = uint16
  uint32_t {.header: "<stdint.h>" importc: "uint32_t".} = uint32
  uint64_t {.header: "<stdint.h>" importc: "uint64_t".} = uint64
  int8_t   {.header: "<stdint.h>" importc: "int8_t".} = int8
  int16_t  {.header: "<stdint.h>" importc: "int16_t".} = int16
  int32_t  {.header: "<stdint.h>" importc: "int32_t".} = int32
  int64_t  {.header: "<stdint.h>" importc: "int64_t".} = int64
type
  abpoa_cigar_t = uint64_t
# const
#   abpoa_cigar_t* = uint64_t

type
  abpoa_res_t* {.bycopy.} = object
    nCigar*: cint
    graphCigar*: ptr abpoa_cigar_t
    nodeS*: cint
    nodeE*: cint
    queryS*: cint
    queryE*: cint             ##  for local and  extension mode
    nAlnBases*: cint
    nMatchedBases*: cint
    bestScore*: int32_t
    isRc* {.bitsize: 1.}: uint8_t
    ##  is_rc: best_score is from the reverse complement

  abpoa_para_t* {.bycopy.} = object
    m*: cint
    mat*: ptr cint              ##  score matrix
    match*: cint
    mismatch*: cint
    gapOpen1*: cint
    gapOpen2*: cint
    gapExt1*: cint
    gapExt2*: cint
    infMin*: cint
    wb*: cint
    wf*: cfloat                ##  extra band width
    zdrop*: cint
    endBonus*: cint           ##  from minimap2
    simdFlag*: cint           ##  available SIMD instruction
                              ##  alignment mode
    retCigar* {.bitsize: 1.}: uint8_t
    revCigar* {.bitsize: 1.}: uint8_t
    outMsa* {.bitsize: 1.}: uint8_t
    outMsaHeader* {.bitsize: 1.}: uint8_t
    outCons* {.bitsize: 1.}: uint8_t
    outGfa* {.bitsize: 1.}: uint8_t
    isDiploid* {.bitsize: 1.}: uint8_t
    useReadIds* {.bitsize: 1.}: uint8_t
    ambStrand* {.bitsize: 1.}: uint8_t
    outPog*: cstring
    alignMode*: cint
    gapMode*: cint
    consAgrm*: cint
    minFreq*: cdouble         ##  for multiploid data
    logTable65536*: array[65536, char]
    bitTable16*: array[65536, char]

  abpoa_node_t* {.bycopy.} = object
    nodeId*: cint
    inEdgeN*: cint
    inEdgeM*: cint
    inId*: ptr cint
    outEdgeN*: cint
    outEdgeM*: cint
    outId*: ptr cint
    outWeight*: ptr cint
    maxOutId*: cint
    readIds*: ptr uint64_t
    readIdsN*: cint          ##  for multiploid
    alignedNodeN*: cint
    alignedNodeM*: cint
    alignedNodeId*: ptr cint  ##  mismatch; aligned node will have same rank
                                ##  int heaviest_weight,
                                ##  heaviest_out_id; // for consensus
    base*: uint8_t             ##  0~m
                 ##  ID, pos ???

  abpoa_graph_t* {.bycopy.} = object
    node*: ptr abpoa_node_t
    nodeN*: cint
    nodeM*: cint
    indexRankM*: cint
    indexToNodeId*: ptr cint
    nodeIdToIndex*: ptr cint
    nodeIdToMaxPosLeft*: ptr cint
    nodeIdToMaxPosRight*: ptr cint
    nodeIdToMaxRemain*: ptr cint
    nodeIdToMsaRank*: ptr cint
    isTopologicalSorted* {.bitsize: 1.}: uint8_t
    isCalledCons* {.bitsize: 1.}: uint8_t
    isSetMsaRank* {.bitsize: 1.}: uint8_t
    calRTime*: cdouble       ##  for evaluation

  abpoa_simd_matrix_t* {.bycopy.} = object
    sMem*: ptr SIMDi
    sMsize*: uint64_t  ##  qp, DP_HE, dp_f OR qp, DP_H, dp_f :
                        ##  based on (qlen, num_of_value, m, node_n)
    dpBeg*: ptr cint
    dpEnd*: ptr cint
    dpBegSn*: ptr cint
    dpEndSn*: ptr cint
    rangM*: cint       ##  if band : based on (node_m)

  abpoa_t* {.bycopy.} = object
    abg*: ptr abpoa_graph_t
    abm*: ptr abpoa_simd_matrix_t


##  init for abpoa parameters

proc abpoa_init_para*(): ptr abpoa_para_t {.importc: "abpoa_init_para".}


proc abpoa_post_set_para*(abpt: ptr abpoa_para_t)
  {.importc: "abpoa_post_set_para".}


proc abpoa_free_para*(abpt: ptr abpoa_para_t)
  {.importc: "abpoa_free_para".}
##  init for alignment


proc abpoa_init*(): ptr abpoa_t
  {.importc: "abpoa_init".}


proc abpoa_free*(ab: ptr abpoa_t;
                 abpt: ptr abpoa_para_t)
                 {.importc: "abpoa_free".}
##  perform msa


proc abpoa_msa*(ab: ptr abpoa_t;
                abpt: ptr abpoa_para_t;
                n_seqs: cint;
                seq_names: cstringArray;
                seq_lens: ptr cint;
                seqs: ptr ptr uint8_t;
                out_fp: ptr File;
                cons_seq: ptr ptr ptr uint8_t;
                cons_cov: ptr ptr ptr cint;
                cons_l: ptr ptr cint;
                cons_n: ptr cint;
                msa_seq: ptr ptr ptr uint8_t;
                msa_l: ptr cint) :
                cint
                {.importc: "abpoa_msa".}
##  clean alignment graph


proc abpoa_reset_graph*(ab: ptr abpoa_t;
                        abpt: ptr abpoa_para_t;
                        qlen: cint)
                        {.importc: "abpoa_reset_graph".}
##  for development:
##  align a sequence to a graph


proc abpoa_align_sequence_to_graph*(ab: ptr abpoa_t;
                                    abpt: ptr abpoa_para_t;
                                    query: ptr uint8_t;
                                    qlen: cint;
                                    res: ptr abpoa_res_t) :
                                    cint
                                    {.importc: "abpoa_align_sequence_to_graph".}
##  align a sequence to a graph between beg_node_id and end_node_id
##    (both are excluded)


proc abpoa_subgraph_nodes*(ab: ptr abpoa_t;
                          inc_beg: cint;
                          inc_end: cint;
                          exc_beg: ptr cint;
                          exc_end: ptr cint)
                          {.importc: "abpoa_subgraph_nodes".}


proc abpoa_align_sequence_to_subgraph*(ab: ptr abpoa_t;
                                       abpt: ptr abpoa_para_t;
                                       beg_node_id: cint;
                                       end_node_id: cint;
                                       query: ptr uint8_t;
                                       qlen: cint;
                                       res: ptr abpoa_res_t) :
                                       cint
                                       {.importc:
                                         "abpoa_align_sequence_to_subgraph".}
##  add a node to a graph
##  para:
##    base: 0123 for ACGT


proc abpoa_add_graph_node*(abg: ptr abpoa_graph_t;
                           base: uint8_t) :
                           cint
                           {.importc: "abpoa_add_graph_node".}
##  add an edge to a graph
##  para:
##    from_id/to_id: ids of from and to nodes
##    check_edge: set as 1 if this edge maybe alread exist and only need to 
##                update weight, set as 0 if the edge is new
##    add_read_id: set as 1 if read_id is used (to use row-column algorithm/
##                 generate MSA result/diploid consensus)
##    read_id: is of sequence
##    read_ids_n: size of read_id array, each one is 64-bit
##                (1+(tot_read_n-1)/64)


proc abpoa_add_graph_edge*(abg: ptr abpoa_graph_t;
                           from_id: cint;
                           to_id: cint;
                           check_edge: cint;
                           w: cint;
                           add_read_id: uint8_t;
                          read_id: cint;
                          read_ids_n: cint) :
                          cint
                          {.importc: "abpoa_add_graph_edge".}
##  add an alignment to a graph
##  para:
##    query: 0123 for ACGT
##    qlen: query length
##    nCigar/abpoa_cigar: from alignment result (abpoa_res_t)
##    read_id: id of sequence
##    tot_read_n: total number of sequence


proc abpoa_add_graph_alignment*(ab: ptr abpoa_t;
                                abpt: ptr abpoa_para_t;
                                query: ptr uint8_t;
                                qlen: cint;
                                res: abpoa_res_t;
                                read_id: cint;
                                tot_read_n: cint) : 
                                cint
                                {.importc: "abpoa_add_graph_alignment".}


proc abpoa_add_subgraph_alignment*(ab: ptr abpoa_t;
                                   abpt: ptr abpoa_para_t;
                                   beg_node_id: cint;
                                   end_node_id: cint;
                                   query: ptr uint8_t;
                                   qlen: cint; res: abpoa_res_t;
                                   read_id: cint;
                                   tot_read_n: cint) : 
                                   cint
                                   {.importc: "abpoa_add_subgraph_alignment".}


proc abpoa_BFS_set_node_index*(abg: ptr abpoa_graph_t;
                               src_id: cint;
                               sink_id: cint) 
                               {.importc: "abpoa_BFS_set_node_index".}


proc abpoa_BFS_set_node_remain*(abg: ptr abpoa_graph_t;
                                src_id: cint;
                                sink_id: cint)
                                {.importc: "abpoa_BFS_set_node_remain".}
##  topological sortting of graph


proc abpoa_topological_sort*(abg: ptr abpoa_graph_t;
                             abpt: ptr abpoa_para_t)
                             {.importc: "abpoa_topological_sort".}
##  generate consensus sequence from graph
##  para:
##    out_fp: consensus sequence output in FASTA format, set as NULL to disable
##    cons_seq, cons_l, cons_n: store consensus sequences in variables,
##                              set cons_n as NULL to disable.
##      cons_seq: store consensus sequences
##      cons_l: store consensus sequences length
##      cons_n: store number of consensus sequences
##      Note: cons_seq and cons_l need to be freed by user.


proc abpoa_generate_consensus*(ab: ptr abpoa_t;
                               abpt: ptr abpoa_para_t;
                               seq_n: cint;
                               out_fp: ptr File;
                               cons_seq: ptr ptr ptr uint8_t;
                               cons_cov: ptr ptr ptr cint;
                               cons_l: ptr ptr cint;
                               cons_n: ptr cint) :
                               cint
                               {.importc: "abpoa_generate_consensus".}
##  generate column multiple sequence alignment from graph


proc abpoa_generate_rc_msa*(ab: ptr abpoa_t;
                            abpt: ptr abpoa_para_t;
                            read_names: cstringArray;
                            is_rc: ptr uint8_t;
                            seq_n: cint;
                            out_fp: ptr File;
                            msa_seq: ptr ptr ptr uint8_t;
                            msa_l: ptr cint)
                            {.importc: "abpoa_generate_rc_msa".}
##  generate graph in GFA format to _out_fp_


proc abpoa_generate_gfa*(ab: ptr abpoa_t;
                         abpt: ptr abpoa_para_t;
                         read_names: cstringArray;
                         is_rc: ptr uint8_t;
                         seq_n: cint;
                         out_fp: ptr File)
                         {.importc: "abpoa_generate_gfa".}
##  generate DOT graph plot and dump graph into PDF/PNG format file


proc abpoa_dump_pog*(ab: ptr abpoa_t;
                     abpt: ptr abpoa_para_t):
                      cint {.importc: "abpoa_dump_pog".}