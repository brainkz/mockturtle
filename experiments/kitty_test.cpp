
#include <kitty/kitty.hpp>
#include <iostream>
#include <vector>
#include <tuple>
#include <fmt/format.h>

constexpr uint8_t NUM_VARS = 4u;

int main()
{
    kitty::static_truth_table<4u> a, b, c, d;
    a._bits = 0x5555;
    b._bits = 0x3333;
    c._bits = 0x0F0F;
    d._bits = 0x00FF;
    std::vector<std::tuple<std::string, kitty::static_truth_table<4u>>> tts = { std::make_tuple("a", a), 
                                                                                std::make_tuple("b", b),
                                                                                std::make_tuple("c", c),
                                                                                std::make_tuple("d", d) } ;
    
    for (uint8_t i = 0u; i < NUM_VARS; ++i)
    {
        for ( auto & [name, tt] : tts )
        {
            if (kitty::has_var(tt, i))
            {
                fmt::print("{0} = {1} = {2} depends on variable #{3}\n", name, kitty::to_hex(tt), kitty::to_binary(tt), i);
            }
        }
    }
}