
add_executable(flimanalysis
        pt3io.c
        utils.c
        lm_mle_diff.c
        fit.c
        convolution.c
        gnuplot_i.c
        em.c
        fileio.c
        readpt3.c
        read_qadata.c
        ini.c
        ferrorf.c
        chgext.c
        fitpt3data.c)
target_link_libraries(flimanalysis GSL::gsl ${CONAN_LIBS})