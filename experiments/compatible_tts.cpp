#include <iostream>
#include <tuple>
#include <vector>
#include <unordered_map>

#include <mockturtle/lib/kitty/kitty.hpp>



constexpr uint8_t NUM_VARS = 3u;

int main()
{
  XOR3._bits = 0x96;
  XNOR3._bits = 0x69;
  MAJ3._bits = 0xe8;
  NMAJ3._bits = 0x17;
  OR3._bits  = 0xfe;
  NOR3._bits  = 0x01;
}
