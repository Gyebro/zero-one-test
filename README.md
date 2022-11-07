# Test for Chaos: The 0-1 test
## Application example using the micro-chaos map

This is a simple C++ implementation of the 0-1 test (See [Georg A. Gottwald and Ian Melbourne; The 0-1 Test for Chaos: A review](https://www.maths.usyd.edu.au/u/gottwald/preprints/testforchaos_MPI.pdf).)
The system subjected to the test is the so-called Micro-Chaos map (see [G. Haller and G. Stépán; Micro-Chaos in Digital Control](https://www.mm.bme.hu/~stepan/docs/Micro_chaos_jnls_1996.pdf).)

The repository uses CMake and relies on GnuPlot.

### Configuring (using CMake)
```
mkdir cmake-build-release
cd cmake-build-release
cmake -DCMAKE_BUILD_TYPE=Release ../ 
cd ..
```
### Building
```
cmake --build cmake-build-release --target zero_one_test -- -j 4
```

### Example output
During the 0-1 test, a sample `Dc-n` plot is generated:
![Example output Dc versus n](/docs/test_for_chaos_snapshot.png)
Then all correlation coefficients (with different `c` values) are plotted on a `Kc-c` graph, like this one:
![Example output Kc versus c](/docs/test_for_chaos_Kc-c.png)
