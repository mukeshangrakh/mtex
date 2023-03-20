import matlab.unittest.plugins.TestReportPlugin;

runner = matlab.unittest.TestRunner.withTextOutput;
runner.addPlugin(TestReportPlugin.producingHTML('Verbosity',3))
runner.run(testsuite({'tConvolution.m', 'tConvolutionFunc.m'}))