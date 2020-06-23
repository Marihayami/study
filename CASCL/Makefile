#Makefile
all:SCFlip_Func.o SCFlip_Main.o
	clang++ -o3 -std=c++11 -o all SCFlip_Func.o SCFlip_Main.o
SCFlip_Func.o:SCFlip_Func.cpp SCFlip.hpp
	clang++ -std=c++11 -c SCFlip_Func.cpp
SCFlip_Main.o:SCFlip_Main.cpp SCFlip.hpp
	clang++ -std=c++11 -c SCFlip_Main.cpp