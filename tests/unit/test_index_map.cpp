#include <gtest/gtest.h>

#include "fg/core/error.hpp"
#include "fg/core/index_map.hpp"

TEST(IndexMap, BijectionAllOrders2D) {
  const fg::IndexMap map(3, 2);
  EXPECT_EQ(map.size(), 7 * 5);
  for (int n = -3; n <= 3; ++n)
    for (int m = -2; m <= 2; ++m) {
      const int r = map.index(n, m);
      const auto [n2, m2] = map.orders(r);
      EXPECT_EQ(n, n2);
      EXPECT_EQ(m, m2);
    }
}

TEST(IndexMap, OneDimensional) {
  const fg::IndexMap map(5, 0);
  EXPECT_FALSE(map.is2D());
  EXPECT_EQ(map.size(), 11);
  for (int n = -5; n <= 5; ++n) EXPECT_EQ(map.index(n, 0), n + 5);
  EXPECT_EQ(map.center(), 5);
}

TEST(IndexMap, OutOfRangeThrows) {
  const fg::IndexMap map(2, 2);
  EXPECT_THROW(map.index(3, 0), fg::Error);
  EXPECT_THROW(map.orders(map.size()), fg::Error);
}
