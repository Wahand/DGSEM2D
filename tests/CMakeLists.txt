ENABLE_TESTING ()

ADD_EXECUTABLE (test_legendreNodesAndWeights test_legendreNodesAndWeights.c)

TARGET_LINK_LIBRARIES (test_legendreNodesAndWeights DGSEM2D-LinAdv)

ADD_TEST (NAME TestLegendreNodesAndWeights COMMAND test_legendreNodesAndWeights)
