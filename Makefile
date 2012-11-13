F77=gfortran

strahlgrid: src/strahlgrid.f
	$(F77) -o $@ $^

clean:
	rm -f strahlgrid
