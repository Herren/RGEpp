CC           = clang++
ODEINTPATH   = odeint/include
EIGENPATH    = eigen
RGEPPPATH    = $(PWD)
FLAVOUROBJ   = src/ckm.o src/pmns.o
FLAGS        = -std=c++11 -I $(ODEINTPATH) -I $(EIGENPATH) -I $(RGEPPPATH)/include -I $(RGEPPPATH)/models

mssmtest: examples/mssmtest.o models/mssm.o $(FLAVOUROBJ)
	$(CC) -o examples/mssmtest models/mssm.o examples/mssmtest.o $(FLAVOUROBJ)

smtest: examples/smtest.o models/sm.o $(FLAVOUROBJ)
	$(CC) -o examples/smtest models/sm.o examples/smtest.o $(FLAVOUROBJ)

nusmtest: examples/nusmtest.o models/nusm.o $(FLAVOUROBJ) examples/rundown.o
	$(CC) -o examples/nusmtest models/nusm.o examples/nusmtest.o examples/rundown.o $(FLAVOUROBJ)

numssmtest: examples/numssmtest.o models/numssm.o $(FLAVOUROBJ) examples/rundown.o
	$(CC) -o examples/numssmtest models/numssm.o examples/numssmtest.o examples/rundown.o $(FLAVOUROBJ)

sm_example: examples/sm_example.o models/sm.o $(FLAVOUROBJ)
	$(CC) -o examples/sm_example models/sm.o examples/sm_example.o $(FLAVOUROBJ)

nsm_example: examples/nsm_example.o models/nsm.o $(FLAVOUROBJ)
	$(CC) -o examples/nsm_example models/nsm.o examples/nsm_example.o $(FLAVOUROBJ)

numssm_example: examples/numssm_example.o models/numssm.o $(FLAVOUROBJ) examples/rundown.o
	$(CC) -o examples/numssm_example examples/numssm_example.o models/numssm.o examples/rundown.o $(FLAVOUROBJ)

thdm_example: examples/thdm_example.o models/thdmi.o models/thdmii.o models/thdmx.o models/thdmy.o $(FLAVOUROBJ)
	$(CC) -o examples/thdm_example examples/thdm_example.o models/thdmi.o models/thdmii.o models/thdmx.o models/thdmy.o $(FLAVOUROBJ)

running_plot: examples/running_plot.o models/sm.o $(FLAVOUROBJ)
	$(CC) -o examples/running_plot examples/running_plot.o models/sm.o $(FLAVOUROBJ)   

%.o: %.cpp
	$(CC) $(FLAGS) -c -O2 -o $@ $<

clean:
	rm src/*.o
	rm models/*.o
	rm examples/*.o








