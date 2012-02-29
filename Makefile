strahlgrid: src/strahlgrid.f
	$(F77) -i $^ -o $@

clean:
	rm -f strahlgrid
