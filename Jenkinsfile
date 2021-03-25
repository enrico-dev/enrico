node {
  withEnv(['OMP_NUM_THREADS=2',
           'OPENMC_CROSS_SECTIONS=$HOME/endfb71_hdf5/cross_sections.xml' ]) {
          sh 'mvn -B verify'
    stage('Hello') {
          echo 'Hello $OPENMC_CROSS_SECTIONS'
     }
  }
}
