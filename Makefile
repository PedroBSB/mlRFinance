CC = g++

SRC = src
TEST = test

.PHONY: clean

clean:
	rm -rf $(SRC)/*.o
	rm -rf $(SRC)/*.dll
	rm -rf $(SRC)/*.a
