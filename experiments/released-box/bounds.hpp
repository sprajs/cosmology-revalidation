#pragma once
#include <array>
#include <cstddef>
// Literal supplied binary64 endpoints from the frozen source contract.
// No center/width arithmetic or prior substitution occurs in this transport.
namespace released_box_transport {
inline constexpr std::array<std::size_t,46> active_original = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,45,46};
inline constexpr std::array<double,46> lower = {
  0x1.ccc9eecbfb15bp+4, // original 0
  0x1.01889a0275254p+5, // original 1
  0x1.fff694467381dp+4, // original 2
  0x1.009e83e425aeep+5, // original 3
  0x1.0c577318fc504p+5, // original 4
  0x1.00741f212d773p+5, // original 5
  0x1.ed886594af4f1p+4, // original 6
  0x1.ef26e978d4fe0p+4, // original 7
  0x1.efe219652bd3cp+4, // original 8
  0x1.ef16872b020c4p+4, // original 9
  0x1.f809d495182aap+4, // original 10
  0x1.fa1c432ca57a8p+4, // original 11
  0x1.f7f0068db8badp+4, // original 12
  0x1.02a4a8c154c98p+5, // original 13
  0x1.fe844d013a92bp+4, // original 14
  0x1.fb7ae147ae148p+4, // original 15
  0x1.fa3bcd35a8588p+4, // original 16
  0x1.017d21ff2e48ep+5, // original 17
  0x1.f14a2339c0ebfp+4, // original 18
  0x1.f1afb7e90ff97p+4, // original 19
  0x1.ef7e28240b780p+4, // original 20
  0x1.ddbc6a7ef9db2p+4, // original 21
  0x1.e5fd21ff2e48fp+4, // original 22
  0x1.f1d3c36113405p+4, // original 23
  0x1.f2bb2fec56d5ep+4, // original 24
  0x1.05e147ae147aep+5, // original 25
  0x1.f7d3c36113405p+4, // original 26
  0x1.e1e69ad42c3cap+4, // original 27
  0x1.002305532617cp+5, // original 28
  0x1.f8cccccccccccp+4, // original 29
  0x1.fabb98c7e2824p+4, // original 30
  0x1.ec3d70a3d70a4p+4, // original 31
  0x1.04afec56d5cfbp+5, // original 32
  0x1.fdd844d013a92p+4, // original 33
  0x1.046a161e4f766p+5, // original 34
  0x1.06930be0ded29p+5, // original 35
  0x1.02e2b6ae7d567p+5, // original 36
  -0x1.e1b6e1a8fce3ep-3, // original 37
  -0x1.84e65bea0ba1fp+2, // original 38
  -0x1.7397a8a8ab2f7p-3, // original 39
  0x1.7ae7d566cf41fp+4, // original 40
  -0x1.5b48ee865af87p-3, // original 41
  -0x1.38212d77318fcp+4, // original 42
  -0x1.4ad5171e29b6bp-1, // original 43
  -0x1.524c688e5ed30p-3, // original 45
  0x1.2178c0053e2d6p+3, // original 46
};
inline constexpr std::array<double,46> upper = {
  0x1.d7fd21ff2e48fp+4, // original 0
  0x1.0cbbcd35a8588p+5, // original 1
  0x1.0cf10cb295e9ep+5, // original 2
  0x1.0a0f27bb2fec6p+5, // original 3
  0x1.1b8aa64c2f838p+5, // original 4
  0x1.07a7525460aa7p+5, // original 5
  0x1.fbeecbfb15b57p+4, // original 6
  0x1.f9b645a1cac08p+4, // original 7
  0x1.fe9a6b50b0f28p+4, // original 8
  0x1.ff6872b020c4ap+4, // original 9
  0x1.0404ea4a8c155p+5, // original 10
  0x1.0d8902de00d1cp+5, // original 11
  0x1.0a875f6fd21ffp+5, // original 12
  0x1.0f487fcb923a2p+5, // original 13
  0x1.071930be0ded3p+5, // original 14
  0x1.049eb851eb852p+5, // original 15
  0x1.0213a92a30553p+5, // original 16
  0x1.0aedc5d638866p+5, // original 17
  0x1.029096bb98c7ep+5, // original 18
  0x1.017bb2fec56d6p+5, // original 19
  0x1.027765fd8adacp+5, // original 20
  0x1.fc22d0e560418p+4, // original 21
  0x1.f463886594af5p+4, // original 22
  0x1.03504816f0069p+5, // original 23
  0x1.0f2027525460bp+5, // original 24
  0x1.0cc28f5c28f5cp+5, // original 25
  0x1.02504816f0069p+5, // original 26
  0x1.ee617c1bda512p+4, // original 27
  0x1.0e896bb98c7e2p+5, // original 28
  0x1.06a3d70a3d70ap+5, // original 29
  0x1.07c432ca57a78p+5, // original 30
  0x1.03e147ae147aep+5, // original 31
  0x1.0fba29c779a6bp+5, // original 32
  0x1.0a484b5dcc63fp+5, // original 33
  0x1.10182a9930be0p+5, // original 34
  0x1.1241205bc01a3p+5, // original 35
  0x1.0a90cb295e9e1p+5, // original 36
  0x1.a367d6a8eea12p-3, // original 37
  -0x1.6ddc1e7967cafp+2, // original 38
  0x1.96a5c7fb2bdadp-3, // original 39
  0x1.90fc504816f01p+4, // original 40
  0x1.0b1d77e00b6dfp-3, // original 41
  -0x1.2f7d566cf41f2p+4, // original 42
  0x1.dee9142b302f5p-3, // original 43
  0x1.1d34c42ceb1aap-5, // original 45
  0x1.33644523f67f4p+3, // original 46
};
}
