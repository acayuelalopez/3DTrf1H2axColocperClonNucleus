import ij.IJ
import ij.ImagePlus
import ij.gui.Roi
import ij.gui.ShapeRoi
import ij.measure.ResultsTable
import ij.plugin.ChannelSplitter
import ij.plugin.Duplicator
import ij.plugin.RGBStackMerge
import org.apache.commons.math3.stat.inference.TTest
import mcib3d.geom.Object3D
import mcib3d.geom.Objects3DPopulation
import mcib3d.image3d.ImageInt

// INPUT UI
//
//#@File(label = "Input Segmentation Files Directory", style = "directory") inputFilesSeg
//#@File(label = " Input Raw Files Directory", style = "directory") inputFilesRaw
//#@File(label = "output directory", style = "directory") outputDir
//#@Integer(label = "Nuclei channel", value = 2) nucleiChannel
//#@Integer(label = "Telomere channel", value = 0) telomereChannel
//#@Integer(label = "Marker channel", value = 1) markerChannel


// IDE
//
//
def inputFilesTrf1 = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2024/2024_8_7_jfiel/output/labels/trf1")
def inputFilesH2ax = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2024/2024_8_7_jfiel/output/labels/h2ax")
def inputFilesNuclei = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2024/2024_8_7_jfiel/output/labels/nuclei")
def inputFilesRaw = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2024/2024_8_7_jfiel/output/images")
def outputDir = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2024/2024_8_7_jfiel/output")
def nucleiChannel = 0.intValue()
def trf1Channel = 2.intValue()
def h2axChannel = 1.intValue()
//def headless = true;
//new ImageJ().setVisible(true);

IJ.log("-Parameters selected: ")
IJ.log("    -Input Seg Files Dir: " + inputFilesTrf1)
IJ.log("    -Input Raw Files Dir: " + inputFilesRaw)
IJ.log("    -output Dir: " + outputDir)
IJ.log("    -Nuclei Channel: " + nucleiChannel)
IJ.log("    -Telomere Channel: " + trf1Channel)
IJ.log("    -Marker Channel: " + h2axChannel)
IJ.log("                                                           ");
/** Get files (images) from input directory */
def listOfFiles = inputFilesRaw.listFiles();
def wtRt = new ResultsTable()
def koRt = new ResultsTable()
def pValueRt = new ResultsTable()
def counterWt = 0.intValue()
def counterKo = 0.intValue()
//WT
def wtColocsN_All = new ArrayList<Double>()
/////Trf1
def wtTrf1N_All = new ArrayList<Double>()
def wtTrf1MeanInt_All = new ArrayList<Double>()
def wtTrf1SumInt_All = new ArrayList<Double>()
def wtTrf1StdInt_All = new ArrayList<Double>()
////H2ax
def wtH2axN_All = new ArrayList<Double>()
def wtH2axMeanInt_All = new ArrayList<Double>()
def wtH2axSumInt_All = new ArrayList<Double>()
def wtH2axStdInt_All = new ArrayList<Double>()
//KO
def koColocsN_All = new ArrayList<Double>()
//////Trf1
def koTrf1N_All = new ArrayList<Double>()
def koTrf1MeanInt_All = new ArrayList<Double>()
def koTrf1SumInt_All = new ArrayList<Double>()
def koTrf1StdInt_All = new ArrayList<Double>()
/////H2ax
def koH2axN_All = new ArrayList<Double>()
def koH2axMeanInt_All = new ArrayList<Double>()
def koH2axSumInt_All = new ArrayList<Double>()
def koH2axStdInt_All = new ArrayList<Double>()


for (def i = 0; i < listOfFiles.length; i++) {
    def counter = 0.intValue()
    def tablePerImage = new ResultsTable();
    /** Create image for each file in the input directory */
    def imps = new ImagePlus(inputFilesRaw.getAbsolutePath() + File.separator + listOfFiles[i].getName())
    def cal = imps.getCalibration()
    IJ.log(imps.getTitle())
    /** Split channels */
    def channels = ChannelSplitter.split(imps)

    /** Get trf1 channel */
    def labelTrf1 = new ImagePlus(inputFilesTrf1.getAbsolutePath() + File.separator + listOfFiles[i].getName().replaceAll(".tif", "_cp_masks.tif"))
    def chTrf1ToMeasure = channels[trf1Channel.intValue()]

    /** Get h2ax channel */
    def labelH2ax = new ImagePlus(inputFilesH2ax.getAbsolutePath() + File.separator + listOfFiles[i].getName().replaceAll(".tif", "_cp_masks.tif"))
    def chH2axToMeasure = channels[h2axChannel.intValue()]

    /** Get nuclei channel */
    def labelNuclei = new ImagePlus(inputFilesNuclei.getAbsolutePath() + File.separator + listOfFiles[i].getName().replaceAll(".tif", "_cp_masks.tif"))
    labelNuclei.setCalibration(cal)
    def chNucleiToMeasure = channels[nucleiChannel.intValue()]
    IJ.log(inputFilesNuclei.getAbsolutePath() + File.separator + listOfFiles[i].getName().replaceAll(".tif", "_cp_masks.tif"))

    // Get Nuclei objects population
    def imgNuclei = ImageInt.wrap(extractCurrentStack(labelNuclei));
    def populationNuclei = new Objects3DPopulation(imgNuclei);
    // Get Nuclei signal
    def signalNuclei = ImageInt.wrap(extractCurrentStack(chNucleiToMeasure));

    // Get Trf1 objects population
    def imgTrf1 = ImageInt.wrap(extractCurrentStack(labelTrf1));
    def populationTrf1 = new Objects3DPopulation(imgTrf1);
    // Get Trf1 signal
    def signalTrf1 = ImageInt.wrap(extractCurrentStack(chTrf1ToMeasure));

    // Get h2ax objects population
    def imgH2ax = ImageInt.wrap(extractCurrentStack(labelH2ax));
    def populationH2ax = new Objects3DPopulation(imgH2ax);
    // Get h2ax signal
    def signalH2ax = ImageInt.wrap(extractCurrentStack(chH2axToMeasure));
    IJ.saveAs(RGBStackMerge.mergeChannels(new ImagePlus[]{labelNuclei, chNucleiToMeasure, labelTrf1, chTrf1ToMeasure, labelH2ax, chH2axToMeasure}, false), "Tiff", outputDir.getAbsolutePath() + File.separator + "merge" + File.separator + listOfFiles[i].getName())



    if (listOfFiles[i].getName().contains("WT")) {
        //WT clon Analysis
        def wtNucleiN = new ArrayList<Double>()
        //TRF1
        def wtTrf1N = new ArrayList<Double>()
        def wtTrf1MeanInt = new ArrayList<Double>()
        def wtTrf1SumInt = new ArrayList<Double>()
        def wtTrf1StdInt = new ArrayList<Double>()

        //H2AX
        def wtH2axN = new ArrayList<Double>()
        def wtH2axMeanInt = new ArrayList<Double>()
        def wtH2axSumInt = new ArrayList<Double>()
        def wtH2axStdInt = new ArrayList<Double>()

        //Colocs trf1-h2ax h2ax-trf1
        def wtColocsN = new ArrayList<Double>()
        wtNucleiN.add(populationNuclei.getNbObjects().doubleValue())

        for (def j = 0.intValue(); j < populationNuclei.getNbObjects(); j++) {
            def counterTrf1 = 0.intValue()
            def counterH2ax = 0.intValue()
            def meanIntTrf1 = new ArrayList<Double>()
            def sumIntTrf1 = new ArrayList<Double>()
            def stdIntTrf1 = new ArrayList<Double>()

            def meanIntH2ax = new ArrayList<Double>()
            def sumIntH2ax = new ArrayList<Double>()
            def stdIntH2ax = new ArrayList<Double>()
            def counterColocs = 0.intValue()

            for (def k = 0.intValue(); k < populationTrf1.getNbObjects(); k++) {
                if (populationNuclei.getObject(j).inside(populationTrf1.getObject(k).getCenterAsPoint())) {
                    counterTrf1++
                    meanIntTrf1.add(populationTrf1.getObject(k).getPixMeanValue(signalTrf1))
                    sumIntTrf1.add(populationTrf1.getObject(k).getIntegratedDensity(signalTrf1))
                    stdIntTrf1.add(populationTrf1.getObject(k).getPixStdDevValue(signalTrf1))
                    wtTrf1MeanInt_All.add(populationTrf1.getObject(k).getPixMeanValue(signalTrf1))
                    wtTrf1StdInt_All.add(populationTrf1.getObject(k).getIntegratedDensity(signalTrf1))
                    wtTrf1SumInt_All.add(populationTrf1.getObject(k).getPixStdDevValue(signalTrf1))
                }
            }
            for (def k = 0.intValue(); k < populationH2ax.getNbObjects(); k++) {
                if (populationNuclei.getObject(j).inside(populationH2ax.getObject(k).getCenterAsPoint())) {
                    counterH2ax++
                    meanIntH2ax.add(populationH2ax.getObject(k).getPixMeanValue(signalH2ax))
                    sumIntH2ax.add(populationH2ax.getObject(k).getIntegratedDensity(signalH2ax))
                    stdIntH2ax.add(populationH2ax.getObject(k).getPixStdDevValue(signalH2ax))
                    wtH2axMeanInt_All.add(populationH2ax.getObject(k).getPixMeanValue(signalH2ax))
                    wtH2axStdInt_All.add(populationH2ax.getObject(k).getIntegratedDensity(signalH2ax))
                    wtH2axSumInt_All.add(populationH2ax.getObject(k).getPixStdDevValue(signalH2ax))
                }
            }
            //Evaluating Colocs
            for (def k = 0.intValue(); k < populationTrf1.getNbObjects(); k++) {
                    for (def l = 0.intValue(); l < populationH2ax.getNbObjects(); l++) {
                        if (populationNuclei.getObject(j).inside(populationTrf1.getObject(k).getCenterAsPoint()) && populationNuclei.getObject(j).inside(populationH2ax.getObject(l).getCenterAsPoint())) {
                            if (populationTrf1.getObject(k).inside(populationH2ax.getObject(l).getCenterAsPoint()))
                                counterColocs++
                    }
                }

            }
            wtColocsN.add(counterColocs.doubleValue())
            wtColocsN_All.add(counterColocs.doubleValue())
            wtTrf1N.add(counterTrf1.doubleValue())
            IJ.log(populationNuclei.getObject(j).getValue() + "-------" + counterTrf1 + "-----colocs:   " + counterColocs)
            wtTrf1N_All.add(counterTrf1.doubleValue())
            wtTrf1MeanInt.add(meanIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
            wtTrf1SumInt.add(sumIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
            wtTrf1StdInt.add(stdIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))

            wtH2axN.add(counterH2ax.doubleValue())
            wtH2axN_All.add(counterH2ax.doubleValue())
            wtH2axMeanInt.add(meanIntH2ax.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
            wtH2axSumInt.add(sumIntH2ax.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
            wtH2axStdInt.add(stdIntH2ax.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
        }
        for (def j = 0.intValue(); j < populationNuclei.getNbObjects(); j++) {
            wtRt.incrementCounter()
            counterWt++
            wtRt.setValue("Image Title", counterWt, listOfFiles[i].getName())
            wtRt.setValue("N of Nuclei per Image", counterWt, populationNuclei.getNbObjects())
            wtRt.setValue("Nucleus Label", counterWt, populationNuclei.getObject(j).getValue())
            wtRt.setValue("Nucleus Volume (microns3)", counterWt, populationNuclei.getObject(j).getVolumeUnit())
            wtRt.setValue("N of TRF1 per Nucleus", counterWt, wtTrf1N.get(j).doubleValue())
            wtRt.setValue("N of H2AX per Nucleus", counterWt, wtH2axN.get(j).doubleValue())
            wtRt.setValue("N of Colocs TRF1-H2Ax", counterWt, wtColocsN.get(j).doubleValue())
            wtRt.setValue("Mean Intensity of TRF1 per Nucleus", counterWt, wtTrf1MeanInt.get(j))
            wtRt.setValue("Mean Intensity of H2AX per Nucleus", counterWt, wtH2axMeanInt.get(j))
            wtRt.setValue("Sum Intensity of TRF1 per Nucleus", counterWt, wtTrf1SumInt.get(j))
            wtRt.setValue("Sum Intensity of H2AX per Nucleus", counterWt, wtH2axSumInt.get(j))
            wtRt.setValue("Std Intensity of TRF1 per Nucleus", counterWt, wtTrf1StdInt.get(j))
            wtRt.setValue("Std Intensity of H2AX per Nucleus", counterWt, wtH2axStdInt.get(j))

        }

    }




    if (listOfFiles[i].getName().contains("KO")) {
        //KO clon Analysis
        def koNucleiN = new ArrayList<Double>()
        //TRF1
        def koTrf1N = new ArrayList<Double>()
        def koTrf1MeanInt = new ArrayList<Double>()
        def koTrf1SumInt = new ArrayList<Double>()
        def koTrf1StdInt = new ArrayList<Double>()

        //H2AX
        def koH2axN = new ArrayList<Double>()
        def koH2axMeanInt = new ArrayList<Double>()
        def koH2axSumInt = new ArrayList<Double>()
        def koH2axStdInt = new ArrayList<Double>()

        //Colocs trf1-h2ax h2ax-trf1
        def koColocsN = new ArrayList<Double>()
        koNucleiN.add(populationNuclei.getNbObjects().doubleValue())

        for (def j = 0.intValue(); j < populationNuclei.getNbObjects(); j++) {
            def counterTrf1 = 0.intValue()
            def counterH2ax = 0.intValue()
            def meanIntTrf1 = new ArrayList<Double>()
            def sumIntTrf1 = new ArrayList<Double>()
            def stdIntTrf1 = new ArrayList<Double>()

            def meanIntH2ax = new ArrayList<Double>()
            def sumIntH2ax = new ArrayList<Double>()
            def stdIntH2ax = new ArrayList<Double>()
            def counterColocs = 0.intValue()

            for (def k = 0.intValue(); k < populationTrf1.getNbObjects(); k++) {
                if (populationNuclei.getObject(j).inside(populationTrf1.getObject(k).getCenterAsPoint())) {
                    counterTrf1++
                    meanIntTrf1.add(populationTrf1.getObject(k).getPixMeanValue(signalTrf1))
                    sumIntTrf1.add(populationTrf1.getObject(k).getIntegratedDensity(signalTrf1))
                    stdIntTrf1.add(populationTrf1.getObject(k).getPixStdDevValue(signalTrf1))
                    koTrf1MeanInt_All.add(populationTrf1.getObject(k).getPixMeanValue(signalTrf1))
                    koTrf1SumInt_All.add(populationTrf1.getObject(k).getIntegratedDensity(signalTrf1))
                    koTrf1StdInt_All.add(populationTrf1.getObject(k).getPixStdDevValue(signalTrf1))
                }
            }
            for (def k = 0.intValue(); k < populationH2ax.getNbObjects(); k++) {
                if (populationNuclei.getObject(j).inside(populationH2ax.getObject(k).getCenterAsPoint())) {
                    counterH2ax++
                    meanIntH2ax.add(populationH2ax.getObject(k).getPixMeanValue(signalH2ax))
                    sumIntH2ax.add(populationH2ax.getObject(k).getIntegratedDensity(signalH2ax))
                    stdIntH2ax.add(populationH2ax.getObject(k).getPixStdDevValue(signalH2ax))
                    koH2axMeanInt_All.add(populationH2ax.getObject(k).getPixMeanValue(signalH2ax))
                    koH2axSumInt_All.add(populationH2ax.getObject(k).getIntegratedDensity(signalH2ax))
                    koH2axStdInt_All.add(populationH2ax.getObject(k).getPixStdDevValue(signalH2ax))
                }
            }

            //Evaluating Colocs
            for (def k = 0.intValue(); k < populationTrf1.getNbObjects(); k++) {
                for (def l = 0.intValue(); l < populationH2ax.getNbObjects(); l++) {
                    if (populationNuclei.getObject(j).inside(populationTrf1.getObject(k).getCenterAsPoint()) && populationNuclei.getObject(j).inside(populationH2ax.getObject(l).getCenterAsPoint())) {
                        if (populationTrf1.getObject(k).inside(populationH2ax.getObject(l).getCenterAsPoint()))
                            counterColocs++
                    }
                }

            }
            koColocsN.add(counterColocs.doubleValue())
            koColocsN_All.add(counterColocs.doubleValue())
            koTrf1N.add(counterTrf1.doubleValue())
            koTrf1N_All.add(counterTrf1.doubleValue())
            koTrf1MeanInt.add(meanIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
            koTrf1SumInt.add(sumIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
            koTrf1StdInt.add(stdIntTrf1.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))

            koH2axN.add(counterH2ax.doubleValue())
            koH2axN_All.add(counterH2ax.doubleValue())
            koH2axMeanInt.add(meanIntH2ax.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
            koH2axSumInt.add(sumIntH2ax.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
            koH2axStdInt.add(stdIntH2ax.stream()
                    .mapToDouble(d -> d)
                    .average()
                    .orElse(0.0))
        }
        for (def j = 0.intValue(); j < populationNuclei.getNbObjects(); j++) {
            koRt.incrementCounter()
            counterKo++
            koRt.setValue("Image Title", counterKo, listOfFiles[i].getName())
            koRt.setValue("Nucleus Label", counterKo, populationNuclei.getObject(j).getValue())
            koRt.setValue("N of Nuclei per Image", counterKo, populationNuclei.getNbObjects())
            koRt.setValue("Nucleus Volume (microns3)", counterKo, populationNuclei.getObject(j).getVolumeUnit())
            koRt.setValue("N of TRF1 per Nucleus", counterKo, koTrf1N.get(j).doubleValue())
            koRt.setValue("N of H2AX per Nucleus", counterKo, koH2axN.get(j).doubleValue())
            koRt.setValue("N of Colocs TRF1-H2Ax", counterKo, koColocsN.get(j).doubleValue())
            koRt.setValue("Mean Intensity of TRF1 per Nucleus", counterKo, koTrf1MeanInt.get(j))
            koRt.setValue("Mean Intensity of H2AX per Nucleus", counterKo, koH2axMeanInt.get(j))
            koRt.setValue("Sum Intensity of TRF1 per Nucleus", counterKo, koTrf1SumInt.get(j))
            koRt.setValue("Sum Intensity of H2AX per Nucleus", counterKo, koH2axSumInt.get(j))
            koRt.setValue("Std Intensity of TRF1 per Nucleus", counterKo, koTrf1StdInt.get(j))
            koRt.setValue("Std Intensity of H2AX per Nucleus", counterKo, koH2axStdInt.get(j))

        }
    }


}
def tTest = new TTest()
IJ.log(wtTrf1N_All.size() + "------" + wtTrf1N_All + "---" + koTrf1N_All.size() + "----" + koTrf1N_All)
pValueRt.setValue("p-value N of TRF1 t-Test", 0, tTest.tTest(wtTrf1N_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koTrf1N_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray()).toString())
pValueRt.setValue("p-value N of Colocs t-Test", 0, tTest.tTest(wtColocsN_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koColocsN_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray()).toString())
pValueRt.setValue("p-value Mean Intensity TRF1 t-Test", 0, tTest.tTest(wtTrf1MeanInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koTrf1MeanInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray()).toString())
pValueRt.setValue("p-value Sum Intensity TRF1 t-Test", 0, tTest.tTest(wtTrf1SumInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koTrf1SumInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray()).toString())
pValueRt.setValue("p-value Std Intensity TRF1 t-Test", 0, tTest.tTest(wtTrf1StdInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koTrf1StdInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray()).toString())
pValueRt.setValue("p-value N of H2AX t-Test", 0, tTest.tTest(wtH2axN_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koH2axN_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray()).toString())
pValueRt.setValue("p-value Mean Intensity H2AX t-Test", 0, tTest.tTest(wtH2axMeanInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koH2axMeanInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray()).toString())
pValueRt.setValue("p-value Sum Intensity H2AX t-Test", 0, tTest.tTest(wtH2axSumInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koH2axSumInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray()).toString())
pValueRt.setValue("p-value Std Intensity H2AX t-Test", 0, tTest.tTest(wtH2axStdInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koH2axStdInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray()).toString())
pValueRt.saveAs(outputDir.getAbsolutePath() + File.separator + "csv" + File.separator + "3D_pvalue_tTest_WT_KO_per_Clon.csv")
wtRt.saveAs(outputDir.getAbsolutePath() + File.separator + "csv" + File.separator + "3DAnalysis_WT_per_Clon.csv")
koRt.saveAs(outputDir.getAbsolutePath() + File.separator + "csv" + File.separator + "3DAnalysis_KO_per_Clon.csv")

ImagePlus extractCurrentStack(ImagePlus plus) {
    // check dimensions
    int[] dims = plus.getDimensions();//XYCZT
    int channel = plus.getChannel();
    int frame = plus.getFrame();
    ImagePlus stack;
    // crop actual frame
    if ((dims[2] > 1) || (dims[4] > 1)) {
        IJ.log("hyperstack found, extracting current channel " + channel + " and frame " + frame);
        def duplicator = new Duplicator();
        stack = duplicator.run(plus, channel, channel, 1, dims[3], frame, frame);
    } else stack = plus.duplicate();

    return stack;
}

static double std(ArrayList<Double> table, double mean) {
    // Step 1:
    double meanDef = mean
    double temp = 0;

    for (int i = 0; i < table.size(); i++) {
        int val = table.get(i);

        // Step 2:
        double squrDiffToMean = Math.pow(val - meanDef, 2);

        // Step 3:
        temp += squrDiffToMean;
    }

    // Step 4:
    double meanOfDiffs = (double) temp / (double) (table.size());

    // Step 5:
    return Math.sqrt(meanOfDiffs);
}