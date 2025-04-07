import requests
from requests.auth import HTTPBasicAuth
from bs4 import BeautifulSoup
import os

# List URLs on the website
url = 'https://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_Hydro/'
reqs = requests.get(url)
soup = BeautifulSoup(reqs.text, 'html.parser')
urls = []
for link in soup.find_all('a'):
    urls.append(link.get('href'))
# Filter to URLs that correspond to hydrologically adjusted elevation
elv_urls = [string for string in urls if "/elv_" in string]
elv_urls = [string.replace("./", url) for string in elv_urls]
#elv_urls = elv_urls[0:2]

# Download the files 
def download_with_auth(url, username, password, filepath):
    """
    Function downloads url using the username and password to a specified filepath
    """
    try:
        response = requests.get(url, auth=HTTPBasicAuth(username, password), stream=True)
        response.raise_for_status()

        with open(filepath, 'wb') as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)
        print(f"Downloaded successfully to {filepath}")
    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")
    except IOError as e:
        print(f"File error: {e}")

if __name__ == "__main__":
    download_directory = "/projects/swot/kmcquil/Global_Hillshade/data/raw/merit/"
    username = "hydrography"
    password = "rivernetwork"
    for url in elv_urls:
        filepath = download_directory + os.path.basename(url)
        download_with_auth(url, username, password, filepath)