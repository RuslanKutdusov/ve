n = 10

lab1:
	g++ --std=c++11 lab1.cpp -o lab1
	./lab1
	cat graphic | gnuplot
	rm data
	rm lab1

lab2: lab2.cpp lab2_graphic
	g++ --std=c++11 lab2.cpp -o lab2
	./lab2 $(n)
	cat lab2_graphic | gnuplot
	rm data
	rm lab2

lab3: lab3.cpp lab3_graphic
	g++ --std=c++11 lab3.cpp -o lab3
	./lab3 $(n)
	cat lab3_graphic | gnuplot
	rm data
	rm lab3
