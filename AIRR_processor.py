import re, sys, getopt, argparse, os
VERBOSE = False


def process_airr_tsv(AIRRFilePath,  outFilePath, metadataPath):
    metadataDict = {}
    metadataLabels = None
    with open(metadataPath, 'r') as metadataFile:
        transcriptList = []
        first_line = True
        for line in metadataFile:
            if first_line:
                metadataLabels = line.rstrip().split(',')
                first_line = False
            else:
                splitTranscriptLine = line.rstrip().split(',')
                transcript, metadataValueList = splitTranscriptLine[0], splitTranscriptLine[1:]
                transcriptList.append(transcript)
                assert len(metadataValueList) == len(metadataLabels)
                metadataDict[transcript] = metadataValueList
        assert metadataLabels != None
                  
        if VERBOSE:
            print('process_AIRR_tsv:', transcriptList)
    
    AIRRFile = open(AIRRFilePath, 'r')
    outFile = open(outFilePath, 'w')
    first_row = True
    contig_id_regex = re.compile(r'-\w+_contig_')
    for line in AIRRFile:
        # Only write the header row to the outFilePath once, during the first sample file read
        if first_row:
            line = '\t'.join([line.rstrip(), '\t'.join(metadataLabels)]) + '\n'
            outFile.write(line)
            first_row = False
        else:
            splitLine = line.rstrip().split('\t', 1)
            
            contig_id = splitLine[0]
            contig_id_split = contig_id_regex.split(contig_id)
            transcript = contig_id_split[0]
            contig_id_contig_num = contig_id_split[1]
            
            remainder = splitLine[1]

            if VERBOSE:
                print('process_AIRR_tsv:', transcript)
            if transcript in transcriptList:
                reconstructedLine = contig_id + '\t' + remainder
                if metadataDict != None:
                    metadataToAdd = '\t'.join(metadataDict[transcript])
                    if VERBOSE:
                        print('process_AIRR_tsv:', metadataToAdd)
                    reconstructedLine = '\t'.join([reconstructedLine, metadataToAdd])
                if VERBOSE:
                    print('process_AIRR_tsv:', reconstructedLine)
                outFile.write(reconstructedLine + '\n')
            else:
                outFile.write(line)
    outFile.close()
    AIRRFile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('AIRR_input', help='path to input AIRR TSV file')
    parser.add_argument('metadata_list', help='CSV with metadata labels, and the metadata to add to each transcript')
    parser.add_argument('-o', '--AIRR_output', help='path to write the AIRR TSV file to', default='AIRR_output.tsv')
    parser.add_argument("-v", "--verbose", help="increase logging verbosity",
                    action="store_true", default=False)
    args = parser.parse_args()
    
    
    
    inputAIRRFilePath = args.AIRR_input
    outputAIRRFilePath = args.AIRR_output
    
    metadataPath = args.metadata_list

    if not os.path.exists(inputAIRRFilePath):
        print('ERROR: the AIRR input TSV does not exist')
        return
    if not os.path.exists(args.metadata_list):
        print('ERROR: the specific metadata list CSV does not exist')
        return
    if args.verbose:
        global VERBOSE
        VERBOSE = True 

    process_airr_tsv(inputAIRRFilePath, outputAIRRFilePath, metadataPath)

if __name__ == "__main__":
        main()