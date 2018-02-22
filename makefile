mcf: main.o io.o fiber.o transport.o header.h
	nvcc main.o io.o fiber.o transport.o -o mcf
%.o: %.cu header.h
	nvcc -c $<
clean:
	rm -rf main.o io.o fiber.o transport.o mcf
