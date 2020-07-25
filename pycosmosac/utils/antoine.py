import re
import numpy as np
import json

def get_antoine_nist(cas_no):
    try:
        from bs4 import BeautifulSoup
        import requests
    except:
        raise ImportError("bs4 and requests needed.")
    url="https://webbook.nist.gov/cgi/cbook.cgi?ID="+cas_no+"&Units=SI&Type=ANTOINE#ANTOINE"
    html_content = requests.get(url).text
    soup = BeautifulSoup(html_content, "lxml")
    #print(soup.prettify()) # print the parsed data of html

    table = soup.find("table", attrs={"class": "data"})
    table_data = table.find_all("tr")  # contains 2 rows

    '''
    # Get all the headings of Lists
    headings = []
    for th in table_data[0].find_all("th"):
        # remove any newlines and extra spaces from left and right
        headings.append(th.text.replace('\n', ' ').strip())
    '''

    T_range = []
    A = []
    B = []
    C = []
    for i in range(1, len(table_data)):
        row_content = []
        for td in table_data[i].find_all("td"):
            row_content.append(td.text)
        tmp = re.findall(r"[-+]?\d*\.\d+|\d+", row_content[0])
        tmp = [float(temp) for temp in tmp]
        T_range.append(tmp)
        A.append(float(re.findall(r"[-+]?\d*\.\d+|\d+",row_content[1])[0]))
        B.append(float(re.findall(r"[-+]?\d*\.\d+|\d+",row_content[2])[0]))
        C.append(float(re.findall(r"[-+]?\d*\.\d+|\d+",row_content[3])[0]))

    '''
    T_range = np.asarray(T_range)
    A = np.asarray(A)
    B = np.asarray(B)
    C = np.asarray(C)
    '''
    return [T_range, A, B, C]


def antoine_to_vapor(data, T):
    T_range = data[0]
    A = data[1]
    B = data[2]
    C = data[3]
    P = []
    for i in range(len(T_range)):
        if T >= T_range[i][0] and T <= T_range[i][1]:
            P.append(10.0**(A[i] - (B[i] / (T + C[i]))) * 100000.0) #Pascal
    if len(P) > 0:
        return(np.average(P))
    else:
        return None


if __name__ == "__main__":
    data = get_antoine_nist("7732-18-5")
    print(antoine_to_vapor(data, 300) - 3554.6837736037555)

    data = [[[250.04, 328.57], [350.14, 466.73], [212.4, 293.02]], [4.022, 4.46988, 4.13377], [1062.64, 1354.913, 1102.878], [-44.93, -5.537, -40.46]]
    P = antoine_to_vapor(data, 298.15)
    print(P - 66910.00086027385)
