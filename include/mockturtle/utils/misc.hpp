
#pragma once

#include <iostream>
#include <string>
#include <set>
#include <fmt/format.h>

std::string repeatString(const std::string& str, int count) 
{
    std::string repeatedStr;
    for (int i = 0; i < count; ++i) 
    {
        repeatedStr += str;
    }
    return repeatedStr;
}

template <typename T1, typename T2>
inline auto generic_max(const T1& a, const T2& b) 
{
  return (a > b) ? a : b;
}

template <typename KeyType, typename ValueType>
void printUnorderedMap(const std::unordered_map<KeyType, ValueType>& map, const int indent_lvl = 0) 
{
  for (const auto& entry : map) 
  {
    fmt::print("{}Key: {} - Value: {}\n", repeatString("\t", indent_lvl), entry.first, entry.second);
  }
}

template <typename ValueType>
void printVector(const std::vector<ValueType>& vec, const int indent_lvl = 0) 
{
  for (const auto& entry : vec) 
  {
    fmt::print("{}{}\n", repeatString("\t", indent_lvl), entry);
  }
}

template <typename elem>
bool haveCommonElements(const std::set<elem>& set1, const std::set<elem>& set2) 
{
    for (const auto& element : set1) 
    {
        if (set2.count(element) > 0) 
        {
            return true;  // Found a common element
        }
    }
    return false;  // No common elements found
}

template <typename Key, typename Value>
Value getMaxValue(const std::unordered_map<Key, Value>& map) 
{
  if (map.empty()) 
  {
    std::cerr << "The map is empty." << std::endl;
    // Return a default-constructed Value or handle the error accordingly
    return Value();
  }
  auto max_entry = std::max_element(map.begin(), map.end(),
    [](const auto& lhs, const auto& rhs) 
    {
        return lhs.second < rhs.second;
    });
  return max_entry->second;
}

template <std::size_t ArrSize>
struct ArrayHash
{
  size_t operator()(const std::array<uint32_t, ArrSize>& arr) const
  {
    size_t seed = 0;
    for (uint32_t value : arr) 
    {
      seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

// std::tuple<int, std::vector<gate>> readGenlibFromFile(const std::string& databasePath) {
//     std::vector<gate> gates;
//     int statusCode = 0;

//     std::ifstream inputFile(databasePath);
//     if (!inputFile.is_open()) {
//         std::cerr << "Error opening input file." << std::endl;
//         statusCode = 1;
//         return std::make_tuple(statusCode, gates);
//     }

//     if (lorina::read_genlib(inputFile, lorina::genlib_reader(gates)) != lorina::return_code::success) {
//         std::cerr << "Error reading genlib." << std::endl;
//         statusCode = 1;
//         return std::make_tuple(statusCode, gates);
//     }

//     return std::make_tuple(statusCode, gates);
// }