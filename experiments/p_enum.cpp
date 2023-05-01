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
  auto [t1, _1, _2] = kitty::exact_p_canonization( tt, [](const kitty::dynamic_truth_table & arg) {fmt::print("P representative: {} = {} = {}\n", kitty::to_binary(arg), kitty::to_hex(arg), arg._bits[0]);});
  fmt::print("The smallest representative is {} = {} = {}\n", kitty::to_binary(t1), kitty::to_hex(t1), t1._bits[0]);
  return 0;
}

#if false

  #include <cstdint>
  #include <iostream>
  #include <unordered_set>

  #include <kitty/kitty.hpp>

  /* compile time constant for the number of variables */
  auto constexpr num_vars = 4;

  int main()
  {
    static_assert( num_vars <= 5, "number of variables is limited to 5" );

    /* truth table type in this example */
    using truth_table = kitty::static_truth_table<num_vars>;

    /* set to store all NPN representatives (dynamic to store bits on heap) */
    kitty::dynamic_truth_table map( truth_table::NumBits );

    /* invert bits: 1 means not classified yet */
    std::transform( map.cbegin(), map.cend(), map.begin(), []( auto word ) { return ~word; } );

    /* hash set to store all NPN classes */
    std::unordered_set<truth_table, kitty::hash<truth_table>> classes;

    /* start from 0 */
    int64_t index = 0;
    truth_table tt;

    while ( index != -1 )
    {
      /* create truth table from index value */
      kitty::create_from_words( tt, &index, &index + 1 );

      /* apply NPN canonization and add resulting representative to set;
        while canonization, mark all encountered truth tables in map
      */
      const auto res = kitty::exact_p_canonization( tt, [&map]( const auto& tt ) { kitty::clear_bit( map, *tt.cbegin() ); } );
      classes.insert( std::get<0>( res ) );

      /* find next non-classified truth table */
      index = find_first_one_bit( map );
    }

    std::cout << "[i] enumerated "
              << map.num_bits() << " functions into "
              << classes.size() << " classes." << std::endl;

    return 0;
  }
#endif