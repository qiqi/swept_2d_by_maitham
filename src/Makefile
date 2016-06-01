CXX = mpicxx
CXX_FLAGS = -O3
CXX_FILES = ClassicDiscretization2D.cpp SweptDiscretization2D.cpp UpPyramid.cpp DownPyramid.cpp VerticalBridge.cpp HorizontalBridge.cpp Main.cpp Socket.cpp ServerSocket.cpp base64.cpp

../bin/streamer:	
	$(CXX) $(CXX_FLAGS) $(CXX_FILES) -o $@
