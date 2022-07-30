import re, sys, getopt, argparse, os
VERBOSE = False

def process_vdj_fasta(fastaFilePath, keptTranscriptsFile, outFile, sampleLabel):
    # build list of transcripts to keep
    transcriptsFile = open(keptTranscriptsFile, 'r')
    transcriptList = []
    for line in transcriptsFile:
        transcriptList.append('>' + line.rstrip().split(',')[0])
    if VERBOSE:
        print('process_vdj_fasta, building transcript list:', transcriptList)
    transcriptsFile.close()

    fastaFile = open(fastaFilePath, 'r')
    keep_transcript = False
    # Given that the sequence_id has the format <transcript>-<sample_id>_contig_<contig_id>,
    # extract the transcript and sample_id
    contig_id_regex = re.compile(r'-(\d+)_')
    for line in fastaFile:
        label_to_use = sampleLabel
        # if the line if the header for a cell
        if '>' in line:
            splitLine = contig_id_regex.split(line.rstrip())
            transcript = splitLine[0]
            # extract the original sample_id inside of the "(\d+)" in the regex statement
            original_label = splitLine[1]
            contig = splitLine[2]
            if VERBOSE:
                print('process_vdj_fasta, found in fasta:', transcript)
            # if the transcript of this cell is in the list of transcripts to keep,
            # write it to the output file with the specified label if given
            if transcript in transcriptList:
                if label_to_use == '':
                    label_to_use = original_label
                reconstructedLine = splitLine[0] + '-' + str(label_to_use) + '_' + contig + '\n'
                outFile.write(reconstructedLine)
                keep_transcript = True
        # if the cell's header has already been written, write its sequence
        elif keep_transcript:
            outFile.write(line)
            keep_transcript = False
    
    fastaFile.close()

def process_vdj_annotation(annotationFilePath, keptTranscriptsFile, outFile, sampleLabel, is_first_sample, metadataLabels = None):
    # if adding new metadata columns, create a dictionary mapping transcripts to their metadata values
    metadataDict = None
    if metadataLabels != None:
        metadataDict = dict()
        
    transcriptsFile = open(keptTranscriptsFile, 'r')
    transcriptList = []
    # for each transcript to keep, add the metadata values for the cell to the metadataDict
    # if metadata columns are being added
    for line in transcriptsFile:
        if metadataDict == None:
            assert len(line.rstrip().split(',')) == 1
            transcriptList.append(line.rstrip())
        else:
            splitTranscriptLine = line.rstrip().split(',')
            transcript, metadataValueList = splitTranscriptLine[0], splitTranscriptLine[1:]
            transcriptList.append(transcript)
            assert len(metadataValueList) == len(metadataLabels)
            metadataDict[transcript] = metadataValueList
    if VERBOSE:
        print('process_vdj_annotation:', transcriptList)
    transcriptsFile.close()
    
    annotationFile = open(annotationFilePath, 'r')
    first_row = True
    # Given that the sequence_id has the format <transcript>-<sample_id>_contig_<contig_id>,
    # extract the transcript
    contig_id_regex = re.compile(r'-\w+_contig_')
    for line in annotationFile:
        label_to_use = sampleLabel
        # Only write the header row to the outfile once, during the first sample file read
        if first_row:
            if is_first_sample:
                if metadataLabels != None:
                    line = ','.join([line.rstrip(), ','.join(metadataLabels)]) + '\n'
                outFile.write(line)
            first_row = False
        else:
            splitLine = line.rstrip().split(',', 4)
            
            # 1st annotation file column: barcode
            transcript_and_sample = splitLine[0].split('-')
            transcript = transcript_and_sample[0]
            original_sample = transcript_and_sample[1]
            
            # 2nd annotation file column: is_cell
            is_cell = splitLine[1]
            
            # 3rd annotation file column: contig_id
            contig_id = splitLine[2]
            contig_id_split = contig_id_regex.split(contig_id)
            contig_id_transcript = contig_id_split[0]
            contig_id_contig_num = contig_id_split[1]

            # 4th annotation file column: high_confidence
            high_confidence = splitLine[3]
            
            # the rest of the columns
            remainder = splitLine[4]

            if VERBOSE:
                print('process_vdj_annotation:', transcript)
            # if keeping the current cell, reconstruct the line using the desired label in the barcode and contig_id columns
            if transcript in transcriptList:
                if label_to_use == '':
                    label_to_use = original_sample
                reconstructedLine = transcript + '-' + label_to_use + ',' + is_cell + ',' + contig_id_transcript + '-' + label_to_use + '_contig_' +contig_id_contig_num + ',' + high_confidence + ',' + remainder
                # add metadata values to the end of the line if applicable
                if metadataDict != None:
                    metadataToAdd = ','.join(metadataDict[transcript])
                    if VERBOSE:
                        print('process_vdj_annotation:', metadataToAdd)
                    reconstructedLine = ','.join([reconstructedLine, metadataToAdd])
                if VERBOSE:
                    print('process_vdj_annotation:', reconstructedLine)
                outFile.write(reconstructedLine + '\n')

    annotationFile.close()

def process_vdj(inRefPaths, metadataList, inFastaPaths, inAnnoPaths, outFastaPath, outAnnotationPath, sampleLabels):
    numSamples = len(sampleLabels)
    # if aggregating annotation CSV files, open the annotation output file and process each input annotation and reference file
    if outAnnotationPath != None:
        outAnnotationFile = open(outAnnotationPath, 'w')
        for i in range(numSamples):
            process_vdj_annotation(inAnnoPaths[i], inRefPaths[i], outAnnotationFile, sampleLabels[i], i == 0, metadataList)
        outAnnotationFile.close()
    # if aggregating FASTA files, open the FASTA output file and process each input FASTA and reference file
    if outFastaPath != None:
        outFastaFile = open(outFastaPath, 'w')
        for i in range(numSamples):
            process_vdj_fasta(inFastaPaths[i], inRefPaths[i], outFastaFile, sampleLabels[i])
        outFastaFile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('config', help='CSV with the input paths and metadata label')
    parser.add_argument('-a', '--annotation_output', help='path to write annotation CSV file to', default='contigs_annotation.csv')
    parser.add_argument('-f', '--fasta_output', help='path to write contigs fasta file to', default='contigs.fasta')
    parser.add_argument("-v", "--verbose", help="increase logging verbosity",
                    action="store_true", default=False)
    args = parser.parse_args()
    
    compareAnnotation = False
    compareFasta = False

    metadataList = None
    
    # Lists of input annotation CSV and fasta files to read, parsed from the consideration file
    inputAnnotationFiles = None
    inputFastaFiles = None
    
    outputFastaFile = None
    outputAnnotationFile = None
    
    sampleLabels = []
    inputReferenceFiles = []
    
    if os.path.exists(args.config):
        # Open the configuration file and parse the first line to check whether
        # to aggregate annotation CSV and/or FASTA files and whether metadata columns
        # should be added
        with open(args.config) as configFile:
            # read first line
            generalConfig = configFile.readline().rstrip().split(',')
            compareAnnotation = generalConfig[0] == 'True'
            compareFasta = generalConfig[1] == 'True'
            if not compareAnnotation and not compareFasta:
                print('ERROR: must aggregate annotations, fastas, or both')
                return
            else:
                if compareAnnotation:
                    inputAnnotationFiles = []
                if compareFasta:
                    inputFastaFiles = []
            
            if len(generalConfig) > 2:
                metadataList = generalConfig[2:]
            # for each vdj run to process, keep track of the sample_id and reference file paths
            # to use, and keep track of the annotation CSV and FASTA paths as appropriate
            for line in configFile:
                splitLine = line.rstrip().split(',')
                numConfigs = len(splitLine)
                
                sampleLabels.append(splitLine[0])
                inputReferenceFiles.append(splitLine[1])
                # if 3 values are given, must aggregate either annotation CSV or FASTA files, not both
                if numConfigs == 3:
                    if compareAnnotation and not compareFasta:
                        inputAnnotationFiles.append(splitLine[2])
                    elif compareFasta and not compareAnnotation:
                        inputFastaFiles.append(splitLine[2])
                    else:
                        print('ERROR: three configuration values given, when four expected')
                        return
                # if 4 values are given, must aggregate both annotation CSV and FASTA files
                elif numConfigs == 4:
                    if compareAnnotation and compareFasta:
                        inputAnnotationFiles.append(splitLine[2])
                        inputFastaFiles.append(splitLine[3])
                    else:
                        print('ERROR: four configuration values given, when three expected')
                        return
                else:
                    print('ERROR: incorrect number of configuration values given', numConfigs)
                    return
    else:
        print('ERROR: the configuration file could not be found')
        return
    if compareAnnotation:
        outputAnnotationFile = args.annotation_output
    if compareFasta:
        outputFastaFile = args.fasta_output 
    if args.verbose:
        global VERBOSE
        VERBOSE = True 

    process_vdj(inputReferenceFiles, metadataList, inputFastaFiles, inputAnnotationFiles, outputFastaFile, outputAnnotationFile, sampleLabels)

if __name__ == "__main__":
        main()