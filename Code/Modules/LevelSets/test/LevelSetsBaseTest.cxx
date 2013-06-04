// --------------------------------------------------------------------------------------
// File:          LevelSetsBaseTest.cxx
// Date:          Jun 4, 2013
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//

#include "gtest/gtest.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST( LevelSetsBaseTest, Basics )
{
	ASSERT_TRUE( true );
}
