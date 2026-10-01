# Build the paper.  The SIAM version needs siamart0516.cls from SIAM's macro distribution
# (CTAN: https://ctan.org/tex-archive/macros/latex/contrib/siam);
# point SIAMDIR at the directory that contains it.
SIAMDIR ?= $(HOME)/texmf/tex/latex/siam
TEX = pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build

all: amsart siam

build:
	mkdir -p build

amsart: build
	$(TEX) main_amsart.tex >/dev/null; $(TEX) main_amsart.tex >/dev/null; $(TEX) main_amsart.tex
	@grep -c "Warning" build/main_amsart.log || true

siam: build
	TEXINPUTS=$(SIAMDIR)//: $(TEX) main_siam.tex >/dev/null; TEXINPUTS=$(SIAMDIR)//: $(TEX) main_siam.tex >/dev/null; TEXINPUTS=$(SIAMDIR)//: $(TEX) main_siam.tex

clean:
	rm -rf build
