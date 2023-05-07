#include <cstdint>
#include <iostream>
#include <unordered_set>

#include <kitty/kitty.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <fmt/format.h>

#pragma GCC diagnostic ignored "-Wunused-parameter"

int next_pow_2(int num)
{
  int power = 0;
  num -= 1;
  while(num > 0)
  {
    power++;
    num >>= 1;
  }
  return power;
}

int main(int argc, char** argv)
{
  std::string in_string = argv[1];
  std::string repr;
  int len = in_string.size();
  // int num_vars = next_pow_2(len);
  int num_vars = next_pow_2(len * 4);
  
  fmt::print("string: {}, len: {}, nextpow: {}\n", in_string, len, num_vars);

  kitty::dynamic_truth_table tt ( num_vars ) ;
  kitty::create_from_hex_string( tt, in_string );
  std::vector<kitty::dynamic_truth_table> seen;
  auto [t1, _1, _2] = kitty::exact_npn_canonization( tt, [&](const kitty::dynamic_truth_table & arg) {
    if (std::find(seen.begin(), seen.end(), arg) == seen.end() )
    {
        fmt::print("NPN representative: {} = {} = {}\n", kitty::to_binary(arg), kitty::to_hex(arg), arg._bits[0]);
        seen.push_back(arg);
    }
    });
  fmt::print("The smallest representative is {} = {} = {}\n", kitty::to_binary(t1), kitty::to_hex(t1), t1._bits[0]);
  return 0;
}