SYSTEM     = x86-64_sles10_4.1
LIBFORMAT  = static_pic

####diretorios
TSPPARSERDIR  = ../tsp-parser
METAHDIR      = ../metaheuristics

CPLEXDIR  = /usr/ilog/cplex
CONCERTDIR = /usr/ilog/concert

#### define o compilador
CPPC = g++ -m64 -O3 -std=c++0x
#############################

#### opcoes de compilacao e includes
CCOPT         = -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -g

TSPINC        = $(TSPPARSERDIR)/include
TSPLIBDIR     = $(TSPPARSERDIR)/lib

METAHINC      = $(METAHDIR)/include
METAHLIBDIR   = $(METAHDIR)/lib

CPLEXINC      = $(CPLEXDIR)/include
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CONCERTINC    = $(CONCERTDIR)/include
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCFLAGS       = $(CCOPT) -I$(CPLEXINC) -I$(CONCERTINC) -I$(TSPINC) -I$(METAHINC)
#############################

#### flags do linker
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex \
	    -L$(CONCERTLIBDIR) -lconcert -lm -pthread \
	    -L$(TSPLIBDIR) -ltsp_parser \
	    -L$(METAHLIBDIR) -lmetah \
#############################

#### diretorios com os source files e com os objs files
SRCDIR = src
OBJDIR = obj
INCDIR = include
#############################

#### lista de todos os srcs e todos os objs
SRCS = $(wildcard $(SRCDIR)/*.cpp)
HDRS = $(wildcard $(INCDIR)/*.h)
OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS))
#############################

#### regra principal, gera o executavel
bnc: $(OBJS) 
	@echo "\033[31m \nLinking all objects files: \033[0m"
	$(CPPC) $(CCFLAGS) $(HDRS) $(OBJS) $(CCLNFLAGS) -o $@ 
############################

#inclui os arquivos de dependencias
-include $(OBJS:.o=.d)

#regra para cada arquivo objeto: compila e gera o arquivo de dependencias do arquivo objeto
#cada arquivo objeto depende do .c e dos headers (informacao dos header esta no arquivo de dependencias gerado pelo compiler)
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo  "\033[31m \nCompiling $<: \033[0m"
	$(CPPC) $(CCFLAGS) -c $< -o $@
	@echo  "\033[32m \ncreating $< dependency file: \033[0m"
	$(CPPC) -MM $(CCFLAGS) $< > $(basename $@).d
	@mv -f $(basename $@).d $(basename $@).d.tmp #proximas tres linhas colocam o diretorio no arquivo de dependencias (g++ nao coloca, surprisingly!)
	@sed -e 's|.*:|$(basename $@).o:|' < $(basename $@).d.tmp > $(basename $@).d
	@rm -f $(basename $@).d.tmp

#delete objetos e arquivos de dependencia
clean:
	@echo "\033[31mcleaning obj directory \033[0m"
	@rm -f $(OBJDIR)/*.o $(OBJDIR)/*.d

rebuild: clean bnc
