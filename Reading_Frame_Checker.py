import pandas as pd


class DelFrame:

    def __init__(self):

        self.refseq_gene_file = (r'Data\refseqgenes.txt')
        self.refseq_df = pd.read_csv(self.refseq_gene_file, sep='\t')

        print(f'RefSeq genes imported to dataframe successfully (Size: {(self.refseq_df.memory_usage().sum()/1000000):2.2f} Megabytes).')

    def state_region(self, event_chr, event_start, event_stop):

        self.event_chr = event_chr
        self.event_start = event_start
        self.event_stop = event_stop

        return f'Stated region: {self.event_chr}:{self.event_start}-{self.event_stop}'

    def ask_region(self):

        print('Enter the chromosome (chrX format and GRCh37 build):')
        event = input()
        self.event_chr = str(event.split(":")[0])
        event_coords = str(event.split(":")[1])
        self.event_start = int(event_coords.split("-")[0])
        self.event_stop = int(event_coords.split("-")[1])

        return f'Input region: {self.event_chr}:{self.event_start}-{self.event_stop}'

    def check_frame_symbols(self, file_name="frame_checks_symbol.txt"):

        self.file_name = file_name

        #New RefSeq DataFrame with only those values on the chromosome of interest
        chr_filter = (self.refseq_df['chrom'] == self.event_chr)
        chr_filt_df = self.refseq_df[chr_filter]

        # Filter to just Genes in which the event is wholly encompassed by Gene Transcript
        event_filter =  ((chr_filt_df['txEnd'] > self.event_start) & 
                        (chr_filt_df['txEnd'] > self.event_stop) &
                        (chr_filt_df['txStart'] < self.event_start) & 
                        (chr_filt_df['txStart'] < self.event_stop))

        event_filt_df = (chr_filt_df[event_filter])

        with open(self.file_name, 'w+') as g:


            for index, transcript in event_filt_df.iterrows():
                exon_starts = transcript['exonStarts'].split(',')
                exon_ends = transcript['exonEnds'].split(',')
                exon_frames = transcript['exonFrames'].split(',')

            
                exon_details = zip(exon_starts,exon_ends, exon_frames)

            
                if transcript['strand'] == '-':
                    exon_details = reversed(tuple(exon_details))
                else:
                    exon_details = tuple(exon_details)


                exon_number = 0
                print(transcript['name2'], '(', transcript['name'],  '):', file=g)


                deletion_started = False
                deletion_ended = False
                exon_before_del = 0
                exon_after_del = 0
                frame_before_del = -2
                frame_after_del = -2

                for exon in exon_details:
                    exon_start, exon_end, exon_frame = exon

                    if exon_frame == "":
                        pass

                    elif int(exon_frame) == -1:
                        exon_number += 1
                        print('#', end='', file=g)
        
                    elif ((int(self.event_start) < int(exon_start)) & (int(self.event_stop) > int(exon_end))):
                        exon_number += 1

                        # if deletion_started == True & deletion_ended == True:
                        #     exon_after_del = exon_number
                        #     frame_after_del = exon_frame
                    
                        deletion_started = True

                        print('*', end='', file=g)
                        
                        
                    else:
                        exon_number += 1

                        if deletion_started == False:
                            exon_before_del = exon_number
                            frame_before_del = exon_frame
                        elif deletion_ended == True:
                            pass
                        else:
                            deletion_ended = True
                            exon_after_del = exon_number
                            frame_after_del = exon_frame

                        if int(exon_frame) == 0:
                            print('|', end='', file=g)
                        elif int(exon_frame) == 1:
                            print('>', end='', file=g)
                        elif int(exon_frame) == 2:
                            print('<', end='', file=g)
                        else:
                            print('@', end='', file=g)

                if (frame_before_del == frame_after_del):
                    print ('\nINFRAME DELETION between exons', exon_before_del, 'and ', exon_after_del, file=g)
                else:
                    print ('\nOUT OF FRAME DELETION between exons', exon_before_del, ' and ', exon_after_del, file=g)

                print('\n', file=g)

        print(f'Written file: {self.file_name}')

        genes_in_event = event_filt_df['name2'].unique().tolist()

        return genes_in_event


    def check_frame_detail(self, file_name="frame_checks_detail.txt"):

        self.file_name = file_name

        #New RefSeq DataFrame with only those values on the chromosome of interest
        chr_filter = (self.refseq_df['chrom'] == self.event_chr)
        chr_filt_df = self.refseq_df[chr_filter]

        # Filter to just Genes in which the event is wholly encompassed by Gene Transcript
        event_filter =  ((chr_filt_df['txEnd'] > self.event_start) & 
                        (chr_filt_df['txEnd'] > self.event_stop) &
                        (chr_filt_df['txStart'] < self.event_start) & 
                        (chr_filt_df['txStart'] < self.event_stop))

        event_filt_df = (chr_filt_df[event_filter])

        with open(self.file_name, 'w+') as f:


            for index, transcript in event_filt_df.iterrows():
                exon_starts = transcript['exonStarts'].split(',')
                exon_ends = transcript['exonEnds'].split(',')
                exon_frames = transcript['exonFrames'].split(',')

                exon_details = zip(exon_starts,exon_ends, exon_frames)

                if transcript['strand'] == '-':
                    exon_details = reversed(tuple(exon_details))
                else:
                    exon_details = tuple(exon_details)

                exon_number = 0
                print(transcript['name2'], '(', transcript['name'],  '):', file=f)

                deletion_started = False
                deletion_ended = False
                exon_before_del = 0
                exon_after_del = 0
                frame_before_del = -2
                frame_after_del = -2

                for exon in exon_details:
                    exon_start, exon_end, exon_frame = exon

                    if exon_frame == "":
                        pass

                    elif int(exon_frame) == -1:
                        exon_number += 1
                        print(exon_number, exon_start, exon_end, exon_frame, ' - Non Coding', file=f)

        
                    elif ((int(self.event_start) < int(exon_start)) & (int(self.event_stop) > int(exon_end))):
                        exon_number += 1

                        # if deletion_started == True & deletion_ended == True:
                        #     exon_after_del = exon_number
                        #     frame_after_del = exon_frame
                    
                        deletion_started = True

                        print(exon_number, exon_start, exon_end, exon_frame, ' - Deleted', file=f)
                        
                        
                    else:
                        exon_number += 1

                        if deletion_started == False:
                            exon_before_del = exon_number
                            frame_before_del = exon_frame
                        elif deletion_ended == True:
                            pass
                        else:
                            deletion_ended = True
                            exon_after_del = exon_number
                            frame_after_del = exon_frame

                        print(exon_number, exon_start, exon_end, exon_frame, file=f)
            
                print('', file=f)

        print(f'Written file: {self.file_name}')

        genes_in_event = event_filt_df['name2'].unique().tolist()

        return genes_in_event


class TestClass:

    def test_one(self):
        d1 = DelFrame()
        d1.state_region('chr2', 54213994, 54380367)
        assert d1.check_frame_detail() == ['ACYP2']

    def test_two(self):
        d2 = DelFrame()
        d2.state_region('chr2', 50400000, 51000000)
        assert d2.check_frame_detail() == ['NRXN1']



if __name__ == "__main__":

    d1 = DelFrame()
    # print(d1.state_region('chr2', 50400000, 51000000))
    print(d1.ask_region())
    d1.check_frame_detail(file_name="frame_checks_detail.txt")
    d1.check_frame_symbols(file_name="frame_checks_symbol.txt")


    

    