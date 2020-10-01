INCLUDES = /usr/local/Cellar/boost/1.73.0/include
HEADERS = $(wildcard *.h)

CXX = g++-10 -std=c++17
DEBUG_FLAGS = -DBOOST_DISABLE_ASSERTS -DNDEBUG 

%.out : %.cc $(HEADERS)
	$(CXX)  $< $(foreach inc,$(INCLUDES),-I $(inc)) -o $@

