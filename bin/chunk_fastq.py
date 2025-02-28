import subprocess
import sys
import os
import tempfile
import timeit
import boto3
import pandas as pd
from urllib.parse import urlparse
S3_BUCKET = "pioneer-scratch"
CHUNK_SIZE = 1000000 #000 # bytes
# When running locally on your desktop machine, set up boto3 to use the AWS SSO credentials
#boto3.setup_default_session(profile_name='general-116981805068')


class S3Url(object):
    # From: https://stackoverflow.com/questions/42641315/s3-urls-get-bucket-name-and-path
    def __init__(self, url):
        self._parsed = urlparse(url, allow_fragments=False)

    @property
    def bucket(self):
        return self._parsed.netloc

    @property
    def key(self):
        if self._parsed.query:
            return self._parsed.path.lstrip("/") + "?" + self._parsed.query
        else:
            return self._parsed.path.lstrip("/")

    @property
    def url(self):
        return self._parsed.geturl()


def process_sample_sheet(sample_sheet_path):

    # Read the sample sheet
    sample_sheet = pd.read_csv(sample_sheet_path)
    # Create a new sample sheet with the same columns
    new_sample_sheet = pd.DataFrame(columns=sample_sheet.columns)
    # Create a new directory to store the downloaded fastq files
    downloaded_s3_files_dir = tempfile.mkdtemp()
    # Download the fastq files from S3
    for index, row in sample_sheet.iterrows():
        # Download the fastq file
        path = row['file']
        if not path.startswith('s3://'):
            print(f"Sequence files in the samples sheet must reside in AWS S3")
            return
        # Get the bucket and key from the s3 path
        s = S3Url(path)
        bucket = s.bucket
        key = s.key # path to the file in the bucket

        s3_dir, fastq_name = key.split('/')[-2:]
        print(f'S3 directory: {s3_dir}, fastq name: {fastq_name}')
        print(f"Downloading fastq file: {key} from S3")
        local_file_name = download_s3_file_to_temp(bucket,key, fastq_name).name
        print(f"Successfully downloaded fastq file: {key} from S3 to {local_file_name}")

        # Chunk the fastq file
        chunk_fastq(local_file_name, f"chunked_files_{index}")
        chunked_files = os.listdir(f"chunked_files_{index}")
        # Load the chunked files to s3, by creating a new S3 folder with the fastq name
        # suffixed with 'chunked'
        s3_chunked_folder_path = fastq_name + '_chunked'
        print(f"Uploading chunked files to S3 in folder: {s3_chunked_folder_path}")
        for chunked_file in chunked_files:
            # Get the path to the chunked file
            path_to_chunked_file = os.path.join(f"chunked_files_{index}", chunked_file)
            upload_file_to_s3(path_to_chunked_file,
                              os.path.join(s3_dir, s3_chunked_folder_path, chunked_file),
                              bucket_name = bucket)
            # Update the new sample sheet with the chunked fastq files
            new_row = row.copy()
            new_row['file'] = os.path.join( 's3://', bucket, s3_dir, s3_chunked_folder_path, chunked_file)
            new_sample_sheet.loc[len(new_sample_sheet)]=new_row


    # Save the new sample sheet
    new_sample_sheet_path =  'new_sample_sheet.csv'
    new_sample_sheet.to_csv(new_sample_sheet_path, index=False)
    return new_sample_sheet_path

def list_files_in_s3_bucket(bucket_name=S3_BUCKET):
    s3 = boto3.client('s3')
    response = s3.list_objects_v2(Bucket=bucket_name)
    print("Files in S3 bucket:")
    for obj in response['Contents']:
        print(obj['Key'])


def upload_file_to_s3(file_name, bucket_file_name, bucket_name=S3_BUCKET):
    s3 = boto3.client('s3')
    s3.upload_file(file_name, bucket_name, bucket_file_name)


def download_s3_file_to_temp(bucket_name, bucket_file_name, temp_file_prefix):
    s3 = boto3.client('s3')

    # Create a temporary file
    temp_file = tempfile.NamedTemporaryFile(prefix=temp_file_prefix, suffix= ".fastq" ,delete=False)

    try:
        # Download the S3 file to the temporary file
        s3.download_fileobj(bucket_name, bucket_file_name, temp_file)
        temp_file.seek(0)  # Go to the beginning of the file
        print(f"File downloaded to temporary location: {temp_file.name}")
    except Exception as e:
        print(f"Error downloading file: {e}")
        temp_file.close()
        raise

    return temp_file


def download_file_from_s3(file_name, local_file_name, bucket_name=S3_BUCKET):
    s3 = boto3.client('s3')
    try:
        s3.download_file(bucket_name, file_name, local_file_name)
    except s3.exceptions.ClientError as e:
        if e.response['Error']['Code'] == '404':
            print(f"File '{file_name}' does not exist")
        return False
    return True

def chunk_fastq(fastq_file, chunked_dir):
    """
    Chunk a fastq file into smaller files of size chunk_size.
    """
    # Count the number of reads in the fastq file
    num_reads = count_reads_using_seqkit(fastq_file)
    # Get the file size
    file_size = os.path.getsize(fastq_file)
    # Convert to gigabytes
    file_size_gb = file_size / CHUNK_SIZE
    # Calculate the average size of a read
    reads_per_gb = int( num_reads / file_size_gb )
    start = timeit.default_timer()
    print(f"Chunking fastq file: {fastq_file} into {chunked_dir} directory")
    seqkit_split_args = ['seqkit', 'split2', fastq_file, '-O', chunked_dir, '-f', '-s', str(reads_per_gb)]
    seqkit_split_process = subprocess.Popen(seqkit_split_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = seqkit_split_process.communicate()
    error_str = error.decode('utf-8')
    if error_str.startswith("Error"):
        print(f"Error: {error_str}")
    else:
        print(f"Successfully chunked fastq file into {chunked_dir} directory")
    end = timeit.default_timer()
    print(f"Time taken to chunk file: {(end-start)/60} minutes to chunk fastq file")

def count_reads_using_seqkit(fastq_file):
    start = timeit.default_timer()
    count_sample_args = ['seqkit','stats', '-T', fastq_file]
    count_process = subprocess.Popen(count_sample_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = count_process.communicate()
    error_str = error.decode('utf-8')
    if error_str.startswith("Error"):
        print(f"Error: {error_str}")
        return
    else:
        print(f"Successfully counted reads in fastq file using seqkit")

    out_string = output.decode('utf-8')
    lines = out_string.split("\n")
    reads = lines[1].split('\t')[3]
    end = timeit.default_timer()
    print(f"Time taken: {(end-start)/60} minutes to count reads in fastq file")
    return int(reads)


if __name__ == "__main__":
    # Given a sample sheet, download the fastq files and chunk them into smaller files.
    # Create a new sample sheet using the new chunked fastq file.
    #fastq_file = sys.argv[1]
    #chunked_dir = sys.argv[2]
    #chunk_fastq(fastq_file, chunked_dir)
    process_sample_sheet(sys.argv[1])

