#include <gtest/gtest.h>

#include <cstddef>

#include "sph_grid.hpp"

// NOLINTNEXTLINE
using namespace halo_pm;

namespace {

class TestSphGrid : public ::testing::Test {};

// NOLINTNEXTLINE
TEST_F(TestSphGrid, xyz) {
  auto [grid_data, pix_num] = GenerateGrid(3);

  float expected_data[768 * 3]{
    0.07207475,  0.07207475,  0.9947916,    //
    -0.07207475, 0.07207475,  0.9947916,    //
    -0.07207475, -0.07207475, 0.9947916,    //
    0.07207475,  -0.07207475, 0.9947916,    //
    0.1876013,   0.07770701,  0.9791666,    //
    0.07770701,  0.1876013,   0.9791666,    //
    -0.07770701, 0.1876013,   0.9791666,    //
    -0.1876013,  0.07770701,  0.9791666,    //
    -0.1876013,  -0.07770701, 0.9791666,    //
    -0.07770701, -0.1876013,  0.9791666,    //
    0.07770701,  -0.1876013,  0.9791666,    //
    0.1876013,   -0.07770701, 0.9791666,    //
    0.2922667,   0.07831264,  0.9531250,    //
    0.2139541,   0.2139541,   0.9531250,    //
    0.07831264,  0.2922667,   0.9531250,    //
    -0.07831264, 0.2922667,   0.9531250,    //
    -0.2139541,  0.2139541,   0.9531250,    //
    -0.2922667,  0.07831264,  0.9531250,    //
    -0.2922667,  -0.07831264, 0.9531250,    //
    -0.2139541,  -0.2139541,  0.9531250,    //
    -0.07831264, -0.2922667,  0.9531250,    //
    0.07831264,  -0.2922667,  0.9531250,    //
    0.2139541,   -0.2139541,  0.9531250,    //
    0.2922667,   -0.07831264, 0.9531250,    //
    0.3919734,   0.07796835,  0.9166666,    //
    0.3322990,   0.2220351,   0.9166666,    //
    0.2220351,   0.3322990,   0.9166666,    //
    0.07796835,  0.3919734,   0.9166666,    //
    -0.07796835, 0.3919734,   0.9166666,    //
    -0.2220351,  0.3322990,   0.9166666,    //
    -0.3322990,  0.2220351,   0.9166666,    //
    -0.3919734,  0.07796835,  0.9166666,    //
    -0.3919734,  -0.07796835, 0.9166666,    //
    -0.3322990,  -0.2220351,  0.9166666,    //
    -0.2220351,  -0.3322990,  0.9166666,    //
    -0.07796835, -0.3919734,  0.9166666,    //
    0.07796835,  -0.3919734,  0.9166666,    //
    0.2220351,   -0.3322990,  0.9166666,    //
    0.3322990,   -0.2220351,  0.9166666,    //
    0.3919734,   -0.07796835, 0.9166666,    //
    0.4873443,   0.07718776,  0.8697916,    //
    0.4396396,   0.2240076,   0.8697916,    //
    0.3489000,   0.3489000,   0.8697916,    //
    0.2240076,   0.4396396,   0.8697916,    //
    0.07718776,  0.4873443,   0.8697916,    //
    -0.07718776, 0.4873443,   0.8697916,    //
    -0.2240076,  0.4396396,   0.8697916,    //
    -0.3489000,  0.3489000,   0.8697916,    //
    -0.4396396,  0.2240076,   0.8697916,    //
    -0.4873443,  0.07718776,  0.8697916,    //
    -0.4873443,  -0.07718776, 0.8697916,    //
    -0.4396396,  -0.2240076,  0.8697916,    //
    -0.3489000,  -0.3489000,  0.8697916,    //
    -0.2240076,  -0.4396396,  0.8697916,    //
    -0.07718776, -0.4873443,  0.8697916,    //
    0.07718776,  -0.4873443,  0.8697916,    //
    0.2240076,   -0.4396396,  0.8697916,    //
    0.3489000,   -0.3489000,  0.8697916,    //
    0.4396396,   -0.2240076,  0.8697916,    //
    0.4873443,   -0.07718776, 0.8697916,    //
    0.5779738,   0.07609170,  0.8125000,    //
    0.5385859,   0.2230895,   0.8125000,    //
    0.4624942,   0.3548842,   0.8125000,    //
    0.3548842,   0.4624942,   0.8125000,    //
    0.2230895,   0.5385859,   0.8125000,    //
    0.07609170,  0.5779738,   0.8125000,    //
    -0.07609170, 0.5779738,   0.8125000,    //
    -0.2230895,  0.5385859,   0.8125000,    //
    -0.3548842,  0.4624942,   0.8125000,    //
    -0.4624942,  0.3548842,   0.8125000,    //
    -0.5385859,  0.2230895,   0.8125000,    //
    -0.5779738,  0.07609170,  0.8125000,    //
    -0.5779738,  -0.07609170, 0.8125000,    //
    -0.5385859,  -0.2230895,  0.8125000,    //
    -0.4624942,  -0.3548842,  0.8125000,    //
    -0.3548842,  -0.4624942,  0.8125000,    //
    -0.2230895,  -0.5385859,  0.8125000,    //
    -0.07609170, -0.5779738,  0.8125000,    //
    0.07609170,  -0.5779738,  0.8125000,    //
    0.2230895,   -0.5385859,  0.8125000,    //
    0.3548842,   -0.4624942,  0.8125000,    //
    0.4624942,   -0.3548842,  0.8125000,    //
    0.5385859,   -0.2230895,  0.8125000,    //
    0.5779738,   -0.07609170, 0.8125000,    //
    0.6631012,   0.07471356,  0.7447916,    //
    0.6298505,   0.2203942,   0.7447916,    //
    0.5650165,   0.3550234,   0.7447916,    //
    0.4718502,   0.4718502,   0.7447916,    //
    0.3550234,   0.5650165,   0.7447916,    //
    0.2203942,   0.6298505,   0.7447916,    //
    0.07471356,  0.6631012,   0.7447916,    //
    -0.07471356, 0.6631012,   0.7447916,    //
    -0.2203942,  0.6298505,   0.7447916,    //
    -0.3550234,  0.5650165,   0.7447916,    //
    -0.4718502,  0.4718502,   0.7447916,    //
    -0.5650165,  0.3550234,   0.7447916,    //
    -0.6298505,  0.2203942,   0.7447916,    //
    -0.6631012,  0.07471356,  0.7447916,    //
    -0.6631012,  -0.07471356, 0.7447916,    //
    -0.6298505,  -0.2203942,  0.7447916,    //
    -0.5650165,  -0.3550234,  0.7447916,    //
    -0.4718502,  -0.4718502,  0.7447916,    //
    -0.3550234,  -0.5650165,  0.7447916,    //
    -0.2203942,  -0.6298505,  0.7447916,    //
    -0.07471356, -0.6631012,  0.7447916,    //
    0.07471356,  -0.6631012,  0.7447916,    //
    0.2203942,   -0.6298505,  0.7447916,    //
    0.3550234,   -0.5650165,  0.7447916,    //
    0.4718502,   -0.4718502,  0.7447916,    //
    0.5650165,   -0.3550234,  0.7447916,    //
    0.6298505,   -0.2203942,  0.7447916,    //
    0.6631012,   -0.07471356, 0.7447916,    //
    0.7417668,   0.07305766,  0.6666666,    //
    0.7132612,   0.2163654,   0.6666666,    //
    0.6573452,   0.3513583,   0.6666666,    //
    0.5761679,   0.4728488,   0.6666666,    //
    0.4728488,   0.5761679,   0.6666666,    //
    0.3513583,   0.6573452,   0.6666666,    //
    0.2163654,   0.7132612,   0.6666666,    //
    0.07305766,  0.7417668,   0.6666666,    //
    -0.07305766, 0.7417668,   0.6666666,    //
    -0.2163654,  0.7132612,   0.6666666,    //
    -0.3513583,  0.6573452,   0.6666666,    //
    -0.4728488,  0.5761679,   0.6666666,    //
    -0.5761679,  0.4728488,   0.6666666,    //
    -0.6573452,  0.3513583,   0.6666666,    //
    -0.7132612,  0.2163654,   0.6666666,    //
    -0.7417668,  0.07305766,  0.6666666,    //
    -0.7417668,  -0.07305766, 0.6666666,    //
    -0.7132612,  -0.2163654,  0.6666666,    //
    -0.6573452,  -0.3513583,  0.6666666,    //
    -0.5761679,  -0.4728488,  0.6666666,    //
    -0.4728488,  -0.5761679,  0.6666666,    //
    -0.3513583,  -0.6573452,  0.6666666,    //
    -0.2163654,  -0.7132612,  0.6666666,    //
    -0.07305766, -0.7417668,  0.6666666,    //
    0.07305766,  -0.7417668,  0.6666666,    //
    0.2163654,   -0.7132612,  0.6666666,    //
    0.3513583,   -0.6573452,  0.6666666,    //
    0.4728488,   -0.5761679,  0.6666666,    //
    0.5761679,   -0.4728488,  0.6666666,    //
    0.6573452,   -0.3513583,  0.6666666,    //
    0.7132612,   -0.2163654,  0.6666666,    //
    0.7417668,   -0.07305766, 0.6666666,    //
    0.8122328,   0.0000000,   0.5833333,    //
    0.7966260,   0.1584587,   0.5833333,    //
    0.7504053,   0.3108280,   0.5833333,    //
    0.6753469,   0.4512524,   0.5833333,    //
    0.5743353,   0.5743353,   0.5833333,    //
    0.4512524,   0.6753469,   0.5833333,    //
    0.3108280,   0.7504053,   0.5833333,    //
    0.1584587,   0.7966260,   0.5833333,    //
    0.0000000,   0.8122328,   0.5833333,    //
    -0.1584587,  0.7966260,   0.5833333,    //
    -0.3108280,  0.7504053,   0.5833333,    //
    -0.4512524,  0.6753469,   0.5833333,    //
    -0.5743353,  0.5743353,   0.5833333,    //
    -0.6753469,  0.4512524,   0.5833333,    //
    -0.7504053,  0.3108280,   0.5833333,    //
    -0.7966260,  0.1584587,   0.5833333,    //
    -0.8122328,  0.0000000,   0.5833333,    //
    -0.7966260,  -0.1584587,  0.5833333,    //
    -0.7504053,  -0.3108280,  0.5833333,    //
    -0.6753469,  -0.4512524,  0.5833333,    //
    -0.5743353,  -0.5743353,  0.5833333,    //
    -0.4512524,  -0.6753469,  0.5833333,    //
    -0.3108280,  -0.7504053,  0.5833333,    //
    -0.1584587,  -0.7966260,  0.5833333,    //
    -0.0000000,  -0.8122328,  0.5833333,    //
    0.1584587,   -0.7966260,  0.5833333,    //
    0.3108280,   -0.7504053,  0.5833333,    //
    0.4512524,   -0.6753469,  0.5833333,    //
    0.5743353,   -0.5743353,  0.5833333,    //
    0.6753469,   -0.4512524,  0.5833333,    //
    0.7504053,   -0.3108280,  0.5833333,    //
    0.7966260,   -0.1584587,  0.5833333,    //
    0.8618552,   0.08488533,  0.5000000,    //
    0.8287346,   0.2513939,   0.5000000,    //
    0.7637662,   0.4082415,   0.5000000,    //
    0.6694466,   0.5494007,   0.5000000,    //
    0.5494007,   0.6694466,   0.5000000,    //
    0.4082415,   0.7637662,   0.5000000,    //
    0.2513939,   0.8287346,   0.5000000,    //
    0.08488533,  0.8618552,   0.5000000,    //
    -0.08488533, 0.8618552,   0.5000000,    //
    -0.2513939,  0.8287346,   0.5000000,    //
    -0.4082415,  0.7637662,   0.5000000,    //
    -0.5494007,  0.6694466,   0.5000000,    //
    -0.6694466,  0.5494007,   0.5000000,    //
    -0.7637662,  0.4082415,   0.5000000,    //
    -0.8287346,  0.2513939,   0.5000000,    //
    -0.8618552,  0.08488533,  0.5000000,    //
    -0.8618552,  -0.08488533, 0.5000000,    //
    -0.8287346,  -0.2513939,  0.5000000,    //
    -0.7637662,  -0.4082415,  0.5000000,    //
    -0.6694466,  -0.5494007,  0.5000000,    //
    -0.5494007,  -0.6694466,  0.5000000,    //
    -0.4082415,  -0.7637662,  0.5000000,    //
    -0.2513939,  -0.8287346,  0.5000000,    //
    -0.08488533, -0.8618552,  0.5000000,    //
    0.08488533,  -0.8618552,  0.5000000,    //
    0.2513939,   -0.8287346,  0.5000000,    //
    0.4082415,   -0.7637662,  0.5000000,    //
    0.5494007,   -0.6694466,  0.5000000,    //
    0.6694466,   -0.5494007,  0.5000000,    //
    0.7637662,   -0.4082415,  0.5000000,    //
    0.8287346,   -0.2513939,  0.5000000,    //
    0.8618552,   -0.08488533, 0.5000000,    //
    0.9090593,   0.0000000,   0.4166666,    //
    0.8915920,   0.1773486,   0.4166666,    //
    0.8398613,   0.3478819,   0.4166666,    //
    0.7558552,   0.5050463,   0.4166666,    //
    0.6428020,   0.6428020,   0.4166666,    //
    0.5050463,   0.7558552,   0.4166666,    //
    0.3478819,   0.8398613,   0.4166666,    //
    0.1773486,   0.8915920,   0.4166666,    //
    0.0000000,   0.9090593,   0.4166666,    //
    -0.1773486,  0.8915920,   0.4166666,    //
    -0.3478819,  0.8398613,   0.4166666,    //
    -0.5050463,  0.7558552,   0.4166666,    //
    -0.6428020,  0.6428020,   0.4166666,    //
    -0.7558552,  0.5050463,   0.4166666,    //
    -0.8398613,  0.3478819,   0.4166666,    //
    -0.8915920,  0.1773486,   0.4166666,    //
    -0.9090593,  0.0000000,   0.4166666,    //
    -0.8915920,  -0.1773486,  0.4166666,    //
    -0.8398613,  -0.3478819,  0.4166666,    //
    -0.7558552,  -0.5050463,  0.4166666,    //
    -0.6428020,  -0.6428020,  0.4166666,    //
    -0.5050463,  -0.7558552,  0.4166666,    //
    -0.3478819,  -0.8398613,  0.4166666,    //
    -0.1773486,  -0.8915920,  0.4166666,    //
    -0.0000000,  -0.9090593,  0.4166666,    //
    0.1773486,   -0.8915920,  0.4166666,    //
    0.3478819,   -0.8398613,  0.4166666,    //
    0.5050463,   -0.7558552,  0.4166666,    //
    0.6428020,   -0.6428020,  0.4166666,    //
    0.7558552,   -0.5050463,  0.4166666,    //
    0.8398613,   -0.3478819,  0.4166666,    //
    0.8915920,   -0.1773486,  0.4166666,    //
    0.9382691,   0.09241144,  0.3333333,    //
    0.9022120,   0.2736830,   0.3333333,    //
    0.8314833,   0.4444371,   0.3333333,    //
    0.7288012,   0.5981117,   0.3333333,    //
    0.5981117,   0.7288012,   0.3333333,    //
    0.4444371,   0.8314833,   0.3333333,    //
    0.2736830,   0.9022120,   0.3333333,    //
    0.09241144,  0.9382691,   0.3333333,    //
    -0.09241144, 0.9382691,   0.3333333,    //
    -0.2736830,  0.9022120,   0.3333333,    //
    -0.4444371,  0.8314833,   0.3333333,    //
    -0.5981117,  0.7288012,   0.3333333,    //
    -0.7288012,  0.5981117,   0.3333333,    //
    -0.8314833,  0.4444371,   0.3333333,    //
    -0.9022120,  0.2736830,   0.3333333,    //
    -0.9382691,  0.09241144,  0.3333333,    //
    -0.9382691,  -0.09241144, 0.3333333,    //
    -0.9022120,  -0.2736830,  0.3333333,    //
    -0.8314833,  -0.4444371,  0.3333333,    //
    -0.7288012,  -0.5981117,  0.3333333,    //
    -0.5981117,  -0.7288012,  0.3333333,    //
    -0.4444371,  -0.8314833,  0.3333333,    //
    -0.2736830,  -0.9022120,  0.3333333,    //
    -0.09241144, -0.9382691,  0.3333333,    //
    0.09241144,  -0.9382691,  0.3333333,    //
    0.2736830,   -0.9022120,  0.3333333,    //
    0.4444371,   -0.8314833,  0.3333333,    //
    0.5981117,   -0.7288012,  0.3333333,    //
    0.7288012,   -0.5981117,  0.3333333,    //
    0.8314833,   -0.4444371,  0.3333333,    //
    0.9022120,   -0.2736830,  0.3333333,    //
    0.9382691,   -0.09241144, 0.3333333,    //
    0.9682458,   0.0000000,   0.2500000,    //
    0.9496412,   0.1888953,   0.2500000,    //
    0.8945425,   0.3705316,   0.2500000,    //
    0.8050669,   0.5379285,   0.2500000,    //
    0.6846531,   0.6846531,   0.2500000,    //
    0.5379285,   0.8050669,   0.2500000,    //
    0.3705316,   0.8945425,   0.2500000,    //
    0.1888953,   0.9496412,   0.2500000,    //
    0.0000000,   0.9682458,   0.2500000,    //
    -0.1888953,  0.9496412,   0.2500000,    //
    -0.3705316,  0.8945425,   0.2500000,    //
    -0.5379285,  0.8050669,   0.2500000,    //
    -0.6846531,  0.6846531,   0.2500000,    //
    -0.8050669,  0.5379285,   0.2500000,    //
    -0.8945425,  0.3705316,   0.2500000,    //
    -0.9496412,  0.1888953,   0.2500000,    //
    -0.9682458,  0.0000000,   0.2500000,    //
    -0.9496412,  -0.1888953,  0.2500000,    //
    -0.8945425,  -0.3705316,  0.2500000,    //
    -0.8050669,  -0.5379285,  0.2500000,    //
    -0.6846531,  -0.6846531,  0.2500000,    //
    -0.5379285,  -0.8050669,  0.2500000,    //
    -0.3705316,  -0.8945425,  0.2500000,    //
    -0.1888953,  -0.9496412,  0.2500000,    //
    -0.0000000,  -0.9682458,  0.2500000,    //
    0.1888953,   -0.9496412,  0.2500000,    //
    0.3705316,   -0.8945425,  0.2500000,    //
    0.5379285,   -0.8050669,  0.2500000,    //
    0.6846531,   -0.6846531,  0.2500000,    //
    0.8050669,   -0.5379285,  0.2500000,    //
    0.8945425,   -0.3705316,  0.2500000,    //
    0.9496412,   -0.1888953,  0.2500000,    //
    0.9812653,   0.09664620,  0.1666666,    //
    0.9435558,   0.2862245,   0.1666666,    //
    0.8695860,   0.4648034,   0.1666666,    //
    0.7621985,   0.6255202,   0.1666666,    //
    0.6255202,   0.7621985,   0.1666666,    //
    0.4648034,   0.8695860,   0.1666666,    //
    0.2862245,   0.9435558,   0.1666666,    //
    0.09664620,  0.9812653,   0.1666666,    //
    -0.09664620, 0.9812653,   0.1666666,    //
    -0.2862245,  0.9435558,   0.1666666,    //
    -0.4648034,  0.8695860,   0.1666666,    //
    -0.6255202,  0.7621985,   0.1666666,    //
    -0.7621985,  0.6255202,   0.1666666,    //
    -0.8695860,  0.4648034,   0.1666666,    //
    -0.9435558,  0.2862245,   0.1666666,    //
    -0.9812653,  0.09664620,  0.1666666,    //
    -0.9812653,  -0.09664620, 0.1666666,    //
    -0.9435558,  -0.2862245,  0.1666666,    //
    -0.8695860,  -0.4648034,  0.1666666,    //
    -0.7621985,  -0.6255202,  0.1666666,    //
    -0.6255202,  -0.7621985,  0.1666666,    //
    -0.4648034,  -0.8695860,  0.1666666,    //
    -0.2862245,  -0.9435558,  0.1666666,    //
    -0.09664620, -0.9812653,  0.1666666,    //
    0.09664620,  -0.9812653,  0.1666666,    //
    0.2862245,   -0.9435558,  0.1666666,    //
    0.4648034,   -0.8695860,  0.1666666,    //
    0.6255202,   -0.7621985,  0.1666666,    //
    0.7621985,   -0.6255202,  0.1666666,    //
    0.8695860,   -0.4648034,  0.1666666,    //
    0.9435558,   -0.2862245,  0.1666666,    //
    0.9812653,   -0.09664620, 0.1666666,    //
    0.9965217,   0.0000000,   0.08333333,   //
    0.9773738,   0.1944117,   0.08333333,   //
    0.9206660,   0.3813523,   0.08333333,   //
    0.8285775,   0.5536378,   0.08333333,   //
    0.7046472,   0.7046472,   0.08333333,   //
    0.5536378,   0.8285775,   0.08333333,   //
    0.3813523,   0.9206660,   0.08333333,   //
    0.1944117,   0.9773738,   0.08333333,   //
    0.0000000,   0.9965217,   0.08333333,   //
    -0.1944117,  0.9773738,   0.08333333,   //
    -0.3813523,  0.9206660,   0.08333333,   //
    -0.5536378,  0.8285775,   0.08333333,   //
    -0.7046472,  0.7046472,   0.08333333,   //
    -0.8285775,  0.5536378,   0.08333333,   //
    -0.9206660,  0.3813523,   0.08333333,   //
    -0.9773738,  0.1944117,   0.08333333,   //
    -0.9965217,  0.0000000,   0.08333333,   //
    -0.9773738,  -0.1944117,  0.08333333,   //
    -0.9206660,  -0.3813523,  0.08333333,   //
    -0.8285775,  -0.5536378,  0.08333333,   //
    -0.7046472,  -0.7046472,  0.08333333,   //
    -0.5536378,  -0.8285775,  0.08333333,   //
    -0.3813523,  -0.9206660,  0.08333333,   //
    -0.1944117,  -0.9773738,  0.08333333,   //
    -0.0000000,  -0.9965217,  0.08333333,   //
    0.1944117,   -0.9773738,  0.08333333,   //
    0.3813523,   -0.9206660,  0.08333333,   //
    0.5536378,   -0.8285775,  0.08333333,   //
    0.7046472,   -0.7046472,  0.08333333,   //
    0.8285775,   -0.5536378,  0.08333333,   //
    0.9206660,   -0.3813523,  0.08333333,   //
    0.9773738,   -0.1944117,  0.08333333,   //
    0.9951847,   0.09801714,  0.0000000,    //
    0.9569403,   0.2902846,   0.0000000,    //
    0.8819212,   0.4713967,   0.0000000,    //
    0.7730104,   0.6343932,   0.0000000,    //
    0.6343932,   0.7730104,   0.0000000,    //
    0.4713967,   0.8819212,   0.0000000,    //
    0.2902846,   0.9569403,   0.0000000,    //
    0.09801714,  0.9951847,   0.0000000,    //
    -0.09801714, 0.9951847,   0.0000000,    //
    -0.2902846,  0.9569403,   0.0000000,    //
    -0.4713967,  0.8819212,   0.0000000,    //
    -0.6343932,  0.7730104,   0.0000000,    //
    -0.7730104,  0.6343932,   0.0000000,    //
    -0.8819212,  0.4713967,   0.0000000,    //
    -0.9569403,  0.2902846,   0.0000000,    //
    -0.9951847,  0.09801714,  0.0000000,    //
    -0.9951847,  -0.09801714, 0.0000000,    //
    -0.9569403,  -0.2902846,  0.0000000,    //
    -0.8819212,  -0.4713967,  0.0000000,    //
    -0.7730104,  -0.6343932,  0.0000000,    //
    -0.6343932,  -0.7730104,  0.0000000,    //
    -0.4713967,  -0.8819212,  0.0000000,    //
    -0.2902846,  -0.9569403,  0.0000000,    //
    -0.09801714, -0.9951847,  0.0000000,    //
    0.09801714,  -0.9951847,  0.0000000,    //
    0.2902846,   -0.9569403,  0.0000000,    //
    0.4713967,   -0.8819212,  0.0000000,    //
    0.6343932,   -0.7730104,  0.0000000,    //
    0.7730104,   -0.6343932,  0.0000000,    //
    0.8819212,   -0.4713967,  0.0000000,    //
    0.9569403,   -0.2902846,  0.0000000,    //
    0.9951847,   -0.09801714, 0.0000000,    //
    0.9965217,   0.0000000,   -0.08333333,  //
    0.9773738,   0.1944117,   -0.08333333,  //
    0.9206660,   0.3813523,   -0.08333333,  //
    0.8285775,   0.5536378,   -0.08333333,  //
    0.7046472,   0.7046472,   -0.08333333,  //
    0.5536378,   0.8285775,   -0.08333333,  //
    0.3813523,   0.9206660,   -0.08333333,  //
    0.1944117,   0.9773738,   -0.08333333,  //
    0.0000000,   0.9965217,   -0.08333333,  //
    -0.1944117,  0.9773738,   -0.08333333,  //
    -0.3813523,  0.9206660,   -0.08333333,  //
    -0.5536378,  0.8285775,   -0.08333333,  //
    -0.7046472,  0.7046472,   -0.08333333,  //
    -0.8285775,  0.5536378,   -0.08333333,  //
    -0.9206660,  0.3813523,   -0.08333333,  //
    -0.9773738,  0.1944117,   -0.08333333,  //
    -0.9965217,  0.0000000,   -0.08333333,  //
    -0.9773738,  -0.1944117,  -0.08333333,  //
    -0.9206660,  -0.3813523,  -0.08333333,  //
    -0.8285775,  -0.5536378,  -0.08333333,  //
    -0.7046472,  -0.7046472,  -0.08333333,  //
    -0.5536378,  -0.8285775,  -0.08333333,  //
    -0.3813523,  -0.9206660,  -0.08333333,  //
    -0.1944117,  -0.9773738,  -0.08333333,  //
    -0.0000000,  -0.9965217,  -0.08333333,  //
    0.1944117,   -0.9773738,  -0.08333333,  //
    0.3813523,   -0.9206660,  -0.08333333,  //
    0.5536378,   -0.8285775,  -0.08333333,  //
    0.7046472,   -0.7046472,  -0.08333333,  //
    0.8285775,   -0.5536378,  -0.08333333,  //
    0.9206660,   -0.3813523,  -0.08333333,  //
    0.9773738,   -0.1944117,  -0.08333333,  //
    0.9812653,   0.09664620,  -0.1666666,   //
    0.9435558,   0.2862245,   -0.1666666,   //
    0.8695860,   0.4648034,   -0.1666666,   //
    0.7621985,   0.6255202,   -0.1666666,   //
    0.6255202,   0.7621985,   -0.1666666,   //
    0.4648034,   0.8695860,   -0.1666666,   //
    0.2862245,   0.9435558,   -0.1666666,   //
    0.09664620,  0.9812653,   -0.1666666,   //
    -0.09664620, 0.9812653,   -0.1666666,   //
    -0.2862245,  0.9435558,   -0.1666666,   //
    -0.4648034,  0.8695860,   -0.1666666,   //
    -0.6255202,  0.7621985,   -0.1666666,   //
    -0.7621985,  0.6255202,   -0.1666666,   //
    -0.8695860,  0.4648034,   -0.1666666,   //
    -0.9435558,  0.2862245,   -0.1666666,   //
    -0.9812653,  0.09664620,  -0.1666666,   //
    -0.9812653,  -0.09664620, -0.1666666,   //
    -0.9435558,  -0.2862245,  -0.1666666,   //
    -0.8695860,  -0.4648034,  -0.1666666,   //
    -0.7621985,  -0.6255202,  -0.1666666,   //
    -0.6255202,  -0.7621985,  -0.1666666,   //
    -0.4648034,  -0.8695860,  -0.1666666,   //
    -0.2862245,  -0.9435558,  -0.1666666,   //
    -0.09664620, -0.9812653,  -0.1666666,   //
    0.09664620,  -0.9812653,  -0.1666666,   //
    0.2862245,   -0.9435558,  -0.1666666,   //
    0.4648034,   -0.8695860,  -0.1666666,   //
    0.6255202,   -0.7621985,  -0.1666666,   //
    0.7621985,   -0.6255202,  -0.1666666,   //
    0.8695860,   -0.4648034,  -0.1666666,   //
    0.9435558,   -0.2862245,  -0.1666666,   //
    0.9812653,   -0.09664620, -0.1666666,   //
    0.9682458,   0.0000000,   -0.2500000,   //
    0.9496412,   0.1888953,   -0.2500000,   //
    0.8945425,   0.3705316,   -0.2500000,   //
    0.8050669,   0.5379285,   -0.2500000,   //
    0.6846531,   0.6846531,   -0.2500000,   //
    0.5379285,   0.8050669,   -0.2500000,   //
    0.3705316,   0.8945425,   -0.2500000,   //
    0.1888953,   0.9496412,   -0.2500000,   //
    0.0000000,   0.9682458,   -0.2500000,   //
    -0.1888953,  0.9496412,   -0.2500000,   //
    -0.3705316,  0.8945425,   -0.2500000,   //
    -0.5379285,  0.8050669,   -0.2500000,   //
    -0.6846531,  0.6846531,   -0.2500000,   //
    -0.8050669,  0.5379285,   -0.2500000,   //
    -0.8945425,  0.3705316,   -0.2500000,   //
    -0.9496412,  0.1888953,   -0.2500000,   //
    -0.9682458,  0.0000000,   -0.2500000,   //
    -0.9496412,  -0.1888953,  -0.2500000,   //
    -0.8945425,  -0.3705316,  -0.2500000,   //
    -0.8050669,  -0.5379285,  -0.2500000,   //
    -0.6846531,  -0.6846531,  -0.2500000,   //
    -0.5379285,  -0.8050669,  -0.2500000,   //
    -0.3705316,  -0.8945425,  -0.2500000,   //
    -0.1888953,  -0.9496412,  -0.2500000,   //
    -0.0000000,  -0.9682458,  -0.2500000,   //
    0.1888953,   -0.9496412,  -0.2500000,   //
    0.3705316,   -0.8945425,  -0.2500000,   //
    0.5379285,   -0.8050669,  -0.2500000,   //
    0.6846531,   -0.6846531,  -0.2500000,   //
    0.8050669,   -0.5379285,  -0.2500000,   //
    0.8945425,   -0.3705316,  -0.2500000,   //
    0.9496412,   -0.1888953,  -0.2500000,   //
    0.9382691,   0.09241144,  -0.3333333,   //
    0.9022120,   0.2736830,   -0.3333333,   //
    0.8314833,   0.4444371,   -0.3333333,   //
    0.7288012,   0.5981117,   -0.3333333,   //
    0.5981117,   0.7288012,   -0.3333333,   //
    0.4444371,   0.8314833,   -0.3333333,   //
    0.2736830,   0.9022120,   -0.3333333,   //
    0.09241144,  0.9382691,   -0.3333333,   //
    -0.09241144, 0.9382691,   -0.3333333,   //
    -0.2736830,  0.9022120,   -0.3333333,   //
    -0.4444371,  0.8314833,   -0.3333333,   //
    -0.5981117,  0.7288012,   -0.3333333,   //
    -0.7288012,  0.5981117,   -0.3333333,   //
    -0.8314833,  0.4444371,   -0.3333333,   //
    -0.9022120,  0.2736830,   -0.3333333,   //
    -0.9382691,  0.09241144,  -0.3333333,   //
    -0.9382691,  -0.09241144, -0.3333333,   //
    -0.9022120,  -0.2736830,  -0.3333333,   //
    -0.8314833,  -0.4444371,  -0.3333333,   //
    -0.7288012,  -0.5981117,  -0.3333333,   //
    -0.5981117,  -0.7288012,  -0.3333333,   //
    -0.4444371,  -0.8314833,  -0.3333333,   //
    -0.2736830,  -0.9022120,  -0.3333333,   //
    -0.09241144, -0.9382691,  -0.3333333,   //
    0.09241144,  -0.9382691,  -0.3333333,   //
    0.2736830,   -0.9022120,  -0.3333333,   //
    0.4444371,   -0.8314833,  -0.3333333,   //
    0.5981117,   -0.7288012,  -0.3333333,   //
    0.7288012,   -0.5981117,  -0.3333333,   //
    0.8314833,   -0.4444371,  -0.3333333,   //
    0.9022120,   -0.2736830,  -0.3333333,   //
    0.9382691,   -0.09241144, -0.3333333,   //
    0.9090593,   0.0000000,   -0.4166666,   //
    0.8915920,   0.1773486,   -0.4166666,   //
    0.8398613,   0.3478819,   -0.4166666,   //
    0.7558552,   0.5050463,   -0.4166666,   //
    0.6428020,   0.6428020,   -0.4166666,   //
    0.5050463,   0.7558552,   -0.4166666,   //
    0.3478819,   0.8398613,   -0.4166666,   //
    0.1773486,   0.8915920,   -0.4166666,   //
    0.0000000,   0.9090593,   -0.4166666,   //
    -0.1773486,  0.8915920,   -0.4166666,   //
    -0.3478819,  0.8398613,   -0.4166666,   //
    -0.5050463,  0.7558552,   -0.4166666,   //
    -0.6428020,  0.6428020,   -0.4166666,   //
    -0.7558552,  0.5050463,   -0.4166666,   //
    -0.8398613,  0.3478819,   -0.4166666,   //
    -0.8915920,  0.1773486,   -0.4166666,   //
    -0.9090593,  0.0000000,   -0.4166666,   //
    -0.8915920,  -0.1773486,  -0.4166666,   //
    -0.8398613,  -0.3478819,  -0.4166666,   //
    -0.7558552,  -0.5050463,  -0.4166666,   //
    -0.6428020,  -0.6428020,  -0.4166666,   //
    -0.5050463,  -0.7558552,  -0.4166666,   //
    -0.3478819,  -0.8398613,  -0.4166666,   //
    -0.1773486,  -0.8915920,  -0.4166666,   //
    -0.0000000,  -0.9090593,  -0.4166666,   //
    0.1773486,   -0.8915920,  -0.4166666,   //
    0.3478819,   -0.8398613,  -0.4166666,   //
    0.5050463,   -0.7558552,  -0.4166666,   //
    0.6428020,   -0.6428020,  -0.4166666,   //
    0.7558552,   -0.5050463,  -0.4166666,   //
    0.8398613,   -0.3478819,  -0.4166666,   //
    0.8915920,   -0.1773486,  -0.4166666,   //
    0.8618552,   0.08488533,  -0.5000000,   //
    0.8287346,   0.2513939,   -0.5000000,   //
    0.7637662,   0.4082415,   -0.5000000,   //
    0.6694466,   0.5494007,   -0.5000000,   //
    0.5494007,   0.6694466,   -0.5000000,   //
    0.4082415,   0.7637662,   -0.5000000,   //
    0.2513939,   0.8287346,   -0.5000000,   //
    0.08488533,  0.8618552,   -0.5000000,   //
    -0.08488533, 0.8618552,   -0.5000000,   //
    -0.2513939,  0.8287346,   -0.5000000,   //
    -0.4082415,  0.7637662,   -0.5000000,   //
    -0.5494007,  0.6694466,   -0.5000000,   //
    -0.6694466,  0.5494007,   -0.5000000,   //
    -0.7637662,  0.4082415,   -0.5000000,   //
    -0.8287346,  0.2513939,   -0.5000000,   //
    -0.8618552,  0.08488533,  -0.5000000,   //
    -0.8618552,  -0.08488533, -0.5000000,   //
    -0.8287346,  -0.2513939,  -0.5000000,   //
    -0.7637662,  -0.4082415,  -0.5000000,   //
    -0.6694466,  -0.5494007,  -0.5000000,   //
    -0.5494007,  -0.6694466,  -0.5000000,   //
    -0.4082415,  -0.7637662,  -0.5000000,   //
    -0.2513939,  -0.8287346,  -0.5000000,   //
    -0.08488533, -0.8618552,  -0.5000000,   //
    0.08488533,  -0.8618552,  -0.5000000,   //
    0.2513939,   -0.8287346,  -0.5000000,   //
    0.4082415,   -0.7637662,  -0.5000000,   //
    0.5494007,   -0.6694466,  -0.5000000,   //
    0.6694466,   -0.5494007,  -0.5000000,   //
    0.7637662,   -0.4082415,  -0.5000000,   //
    0.8287346,   -0.2513939,  -0.5000000,   //
    0.8618552,   -0.08488533, -0.5000000,   //
    0.8122328,   0.0000000,   -0.5833333,   //
    0.7966260,   0.1584587,   -0.5833333,   //
    0.7504053,   0.3108280,   -0.5833333,   //
    0.6753469,   0.4512524,   -0.5833333,   //
    0.5743353,   0.5743353,   -0.5833333,   //
    0.4512524,   0.6753469,   -0.5833333,   //
    0.3108280,   0.7504053,   -0.5833333,   //
    0.1584587,   0.7966260,   -0.5833333,   //
    0.0000000,   0.8122328,   -0.5833333,   //
    -0.1584587,  0.7966260,   -0.5833333,   //
    -0.3108280,  0.7504053,   -0.5833333,   //
    -0.4512524,  0.6753469,   -0.5833333,   //
    -0.5743353,  0.5743353,   -0.5833333,   //
    -0.6753469,  0.4512524,   -0.5833333,   //
    -0.7504053,  0.3108280,   -0.5833333,   //
    -0.7966260,  0.1584587,   -0.5833333,   //
    -0.8122328,  0.0000000,   -0.5833333,   //
    -0.7966260,  -0.1584587,  -0.5833333,   //
    -0.7504053,  -0.3108280,  -0.5833333,   //
    -0.6753469,  -0.4512524,  -0.5833333,   //
    -0.5743353,  -0.5743353,  -0.5833333,   //
    -0.4512524,  -0.6753469,  -0.5833333,   //
    -0.3108280,  -0.7504053,  -0.5833333,   //
    -0.1584587,  -0.7966260,  -0.5833333,   //
    -0.0000000,  -0.8122328,  -0.5833333,   //
    0.1584587,   -0.7966260,  -0.5833333,   //
    0.3108280,   -0.7504053,  -0.5833333,   //
    0.4512524,   -0.6753469,  -0.5833333,   //
    0.5743353,   -0.5743353,  -0.5833333,   //
    0.6753469,   -0.4512524,  -0.5833333,   //
    0.7504053,   -0.3108280,  -0.5833333,   //
    0.7966260,   -0.1584587,  -0.5833333,   //
    0.7417668,   0.07305766,  -0.6666666,   //
    0.7132612,   0.2163654,   -0.6666666,   //
    0.6573452,   0.3513583,   -0.6666666,   //
    0.5761679,   0.4728488,   -0.6666666,   //
    0.4728488,   0.5761679,   -0.6666666,   //
    0.3513583,   0.6573452,   -0.6666666,   //
    0.2163654,   0.7132612,   -0.6666666,   //
    0.07305766,  0.7417668,   -0.6666666,   //
    -0.07305766, 0.7417668,   -0.6666666,   //
    -0.2163654,  0.7132612,   -0.6666666,   //
    -0.3513583,  0.6573452,   -0.6666666,   //
    -0.4728488,  0.5761679,   -0.6666666,   //
    -0.5761679,  0.4728488,   -0.6666666,   //
    -0.6573452,  0.3513583,   -0.6666666,   //
    -0.7132612,  0.2163654,   -0.6666666,   //
    -0.7417668,  0.07305766,  -0.6666666,   //
    -0.7417668,  -0.07305766, -0.6666666,   //
    -0.7132612,  -0.2163654,  -0.6666666,   //
    -0.6573452,  -0.3513583,  -0.6666666,   //
    -0.5761679,  -0.4728488,  -0.6666666,   //
    -0.4728488,  -0.5761679,  -0.6666666,   //
    -0.3513583,  -0.6573452,  -0.6666666,   //
    -0.2163654,  -0.7132612,  -0.6666666,   //
    -0.07305766, -0.7417668,  -0.6666666,   //
    0.07305766,  -0.7417668,  -0.6666666,   //
    0.2163654,   -0.7132612,  -0.6666666,   //
    0.3513583,   -0.6573452,  -0.6666666,   //
    0.4728488,   -0.5761679,  -0.6666666,   //
    0.5761679,   -0.4728488,  -0.6666666,   //
    0.6573452,   -0.3513583,  -0.6666666,   //
    0.7132612,   -0.2163654,  -0.6666666,   //
    0.7417668,   -0.07305766, -0.6666666,   //
    0.6631012,   0.07471356,  -0.7447916,   //
    0.6298505,   0.2203942,   -0.7447916,   //
    0.5650165,   0.3550234,   -0.7447916,   //
    0.4718502,   0.4718502,   -0.7447916,   //
    0.3550234,   0.5650165,   -0.7447916,   //
    0.2203942,   0.6298505,   -0.7447916,   //
    0.07471356,  0.6631012,   -0.7447916,   //
    -0.07471356, 0.6631012,   -0.7447916,   //
    -0.2203942,  0.6298505,   -0.7447916,   //
    -0.3550234,  0.5650165,   -0.7447916,   //
    -0.4718502,  0.4718502,   -0.7447916,   //
    -0.5650165,  0.3550234,   -0.7447916,   //
    -0.6298505,  0.2203942,   -0.7447916,   //
    -0.6631012,  0.07471356,  -0.7447916,   //
    -0.6631012,  -0.07471356, -0.7447916,   //
    -0.6298505,  -0.2203942,  -0.7447916,   //
    -0.5650165,  -0.3550234,  -0.7447916,   //
    -0.4718502,  -0.4718502,  -0.7447916,   //
    -0.3550234,  -0.5650165,  -0.7447916,   //
    -0.2203942,  -0.6298505,  -0.7447916,   //
    -0.07471356, -0.6631012,  -0.7447916,   //
    0.07471356,  -0.6631012,  -0.7447916,   //
    0.2203942,   -0.6298505,  -0.7447916,   //
    0.3550234,   -0.5650165,  -0.7447916,   //
    0.4718502,   -0.4718502,  -0.7447916,   //
    0.5650165,   -0.3550234,  -0.7447916,   //
    0.6298505,   -0.2203942,  -0.7447916,   //
    0.6631012,   -0.07471356, -0.7447916,   //
    0.5779738,   0.07609170,  -0.8125000,   //
    0.5385859,   0.2230895,   -0.8125000,   //
    0.4624942,   0.3548842,   -0.8125000,   //
    0.3548842,   0.4624942,   -0.8125000,   //
    0.2230895,   0.5385859,   -0.8125000,   //
    0.07609170,  0.5779738,   -0.8125000,   //
    -0.07609170, 0.5779738,   -0.8125000,   //
    -0.2230895,  0.5385859,   -0.8125000,   //
    -0.3548842,  0.4624942,   -0.8125000,   //
    -0.4624942,  0.3548842,   -0.8125000,   //
    -0.5385859,  0.2230895,   -0.8125000,   //
    -0.5779738,  0.07609170,  -0.8125000,   //
    -0.5779738,  -0.07609170, -0.8125000,   //
    -0.5385859,  -0.2230895,  -0.8125000,   //
    -0.4624942,  -0.3548842,  -0.8125000,   //
    -0.3548842,  -0.4624942,  -0.8125000,   //
    -0.2230895,  -0.5385859,  -0.8125000,   //
    -0.07609170, -0.5779738,  -0.8125000,   //
    0.07609170,  -0.5779738,  -0.8125000,   //
    0.2230895,   -0.5385859,  -0.8125000,   //
    0.3548842,   -0.4624942,  -0.8125000,   //
    0.4624942,   -0.3548842,  -0.8125000,   //
    0.5385859,   -0.2230895,  -0.8125000,   //
    0.5779738,   -0.07609170, -0.8125000,   //
    0.4873443,   0.07718776,  -0.8697916,   //
    0.4396396,   0.2240076,   -0.8697916,   //
    0.3489000,   0.3489000,   -0.8697916,   //
    0.2240076,   0.4396396,   -0.8697916,   //
    0.07718776,  0.4873443,   -0.8697916,   //
    -0.07718776, 0.4873443,   -0.8697916,   //
    -0.2240076,  0.4396396,   -0.8697916,   //
    -0.3489000,  0.3489000,   -0.8697916,   //
    -0.4396396,  0.2240076,   -0.8697916,   //
    -0.4873443,  0.07718776,  -0.8697916,   //
    -0.4873443,  -0.07718776, -0.8697916,   //
    -0.4396396,  -0.2240076,  -0.8697916,   //
    -0.3489000,  -0.3489000,  -0.8697916,   //
    -0.2240076,  -0.4396396,  -0.8697916,   //
    -0.07718776, -0.4873443,  -0.8697916,   //
    0.07718776,  -0.4873443,  -0.8697916,   //
    0.2240076,   -0.4396396,  -0.8697916,   //
    0.3489000,   -0.3489000,  -0.8697916,   //
    0.4396396,   -0.2240076,  -0.8697916,   //
    0.4873443,   -0.07718776, -0.8697916,   //
    0.3919734,   0.07796835,  -0.9166666,   //
    0.3322990,   0.2220351,   -0.9166666,   //
    0.2220351,   0.3322990,   -0.9166666,   //
    0.07796835,  0.3919734,   -0.9166666,   //
    -0.07796835, 0.3919734,   -0.9166666,   //
    -0.2220351,  0.3322990,   -0.9166666,   //
    -0.3322990,  0.2220351,   -0.9166666,   //
    -0.3919734,  0.07796835,  -0.9166666,   //
    -0.3919734,  -0.07796835, -0.9166666,   //
    -0.3322990,  -0.2220351,  -0.9166666,   //
    -0.2220351,  -0.3322990,  -0.9166666,   //
    -0.07796835, -0.3919734,  -0.9166666,   //
    0.07796835,  -0.3919734,  -0.9166666,   //
    0.2220351,   -0.3322990,  -0.9166666,   //
    0.3322990,   -0.2220351,  -0.9166666,   //
    0.3919734,   -0.07796835, -0.9166666,   //
    0.2922667,   0.07831264,  -0.9531250,   //
    0.2139541,   0.2139541,   -0.9531250,   //
    0.07831264,  0.2922667,   -0.9531250,   //
    -0.07831264, 0.2922667,   -0.9531250,   //
    -0.2139541,  0.2139541,   -0.9531250,   //
    -0.2922667,  0.07831264,  -0.9531250,   //
    -0.2922667,  -0.07831264, -0.9531250,   //
    -0.2139541,  -0.2139541,  -0.9531250,   //
    -0.07831264, -0.2922667,  -0.9531250,   //
    0.07831264,  -0.2922667,  -0.9531250,   //
    0.2139541,   -0.2139541,  -0.9531250,   //
    0.2922667,   -0.07831264, -0.9531250,   //
    0.1876013,   0.07770701,  -0.9791666,   //
    0.07770701,  0.1876013,   -0.9791666,   //
    -0.07770701, 0.1876013,   -0.9791666,   //
    -0.1876013,  0.07770701,  -0.9791666,   //
    -0.1876013,  -0.07770701, -0.9791666,   //
    -0.07770701, -0.1876013,  -0.9791666,   //
    0.07770701,  -0.1876013,  -0.9791666,   //
    0.1876013,   -0.07770701, -0.9791666,   //
    0.07207475,  0.07207475,  -0.9947916,   //
    -0.07207475, 0.07207475,  -0.9947916,   //
    -0.07207475, -0.07207475, -0.9947916,   //
    0.07207475,  -0.07207475, -0.9947916,   //
  };

  ASSERT_EQ(pix_num, 768);
  for (size_t i = 0; i < pix_num * 3; i++) {
    EXPECT_NEAR(expected_data[i], grid_data[i], 1e-6) << "(" << i << ")";
  }
}

}  // namespace
