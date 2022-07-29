import datetime
import os
import tempfile
import hashlib
import pytest
from calicolipidlibrary import *


def file_checksum(fname, chunksize=8192):
    with open(fname, "rb") as f:
        file_hash = hashlib.md5()
        chunk = f.read(chunksize)
        while chunk:
            file_hash.update(chunk)
            chunk = f.read(chunksize)
    return file_hash.hexdigest()


def checksum_for_generated_library(lipid_cls, ion):
    f = tempfile.NamedTemporaryFile(mode='w+', delete=False)
    f.close()
    lipid_cls().generateLibrary(f.name, mode=ion)
    csum = file_checksum(f.name)
    os.unlink(f.name)  # delete the file
    return csum

test_params = [(MIP2C,         "143edab66a6f3f9d05acdf3b411d72bc", "8f23c164cf2bc31279e8dbcd31e85077"),
               ]


test_params2 = [(AcGM2,         "ea411225fdcfc944932266c142443710", "0f72e21f828f93b0ac1bb9ceba9865f2"),
               (AcGM3,         "8178e54d2ee356d7c17a1eacda122835", "57f7bd38ec39732681216b2aebdbbccc"),
               (Alkyl_LPC,     "38fa54a34da19d9505e041379056a9b6", "4edda981da32bd3950ab7bc0bd7f2864"),
               (Alkyl_LPE,     "a576f9fe68ad2abd2f6e8621a39e62a7", "a07f86b35d7679ecc2482b2cc8c0b5f9"),
               (Alkyl_LPS,     "c96a3a1a59dbfe36ae88b4a87895af53", "7e8562de40b11cf4a3e1bdb1a8144eee"),
               (Alkyl_PC,      "4bc6a0b3a5dcdb025b81fad96db4c501", "b90f66dd169c469b8930ca1bfd761c82"),
               (Alkyl_PE,      "a84b83dafed85324a7dced10f0b3b667", "c7085349a6398493c6cf5332b83a179b"),
               (Alkyl_PS,      "22aab1538df291b1e0296a0b6fc2181c", "af897221d1952d72be5651fbf1fa7457"),
               (BDP,           "f485fe58e198388d167cd0c9aecea5dd", "6d4358098013f6665bb90c5e99e12076"),
               (BMP,           "2a17651f16575da65f4b5c3dac0a6ac4", "1c1234f91884035c24705e4ff821a477"),
               (Carn,          "c7814cd5d1f354354a8047f6daf8550f", "bd144d4e286056405e5488e3ad308e24"),
               (CDP_DG,        "093034431b686789540ae148ec53ea63", "15056114de0593ebe6c2927468d6c76b"),
               (CE,            "a7ff6d809525245919d3a56933312ae1", "ee55391fc7c2e68386f1650e115dadb1"),
               (Ceramide,      "63b70f3a0ca4f4f68c463b5d64aed3b7", "15b1af6f54843a3e55693066f9e16bc6"),
               (Ceramide_P,    "5def6a6836240a65a86d029e2c6f9430", "6d10b04cd3edcacbae0f8ccadf17c654"),
               (CL,            "ffeed408f8fbe22fad1e5fef70eaa91e", "f834b408ad54d26c36d023eeaa463a50"),
               (CPE,           "88c8edbeeaf9ec7657b61c8fa1c7281f", "8a04775cc3dca9079041f1c8d291a333"),
               (CPI,           "db1fcb5aca6353a1f457b09b19b52d7b", "feb6d081a616985236fa0433b89afce9"),
               (DG,            "6a81ed4ab7dc163dd7610b06c0193c7e", "d83067139ea3593ec0a1d0111d73e123"),
               (DGDG,          "bdcb42ff2238c4194ffa244f4dd4a1ad", "1297975670fcd889794aae547510613f"),
               (DGTS,          "e29bb15af39c4a3b01a473cacc9ac31c", "a12a016020b8184a543979c5c2b26dbf"),
               (DMPE,          "d837a1895deed1d7909a486e10b799dd", "5572479e1a595c4f8c1356dd77d21209"),
               (ErgE,          "7b5d73e23a179642d02fa042ed39d26b", "36df1d7cff78abc4009fac39b80cbf29"),
               (Ethanolamine,  "eb175c24bcf1bfd3bce128f99ad538bb", "9aefcb594a46cee0fc90fbd998addcb2"),
               (FA,            "39f0a2c82c477ec8bc978bf171befa87", "3246d63d78c01123d3f42835f2d5b85f"),
               (FAHFA,         "242f86493d08178edd63581a49c8407c", "6136a4143144813c4f7e8084e2a46df4"),
               (GB3,           "46d18521bdd7f85a50d5044f405b7fe3", "ff28732f1c77f027427d450288d7b643"),
               (GcGM2,         "99c8303812627ca3f63c12b43efd5cd9", "81910925ea8257b8ef029f7e7bab5abe"),
               (GcGM3,         "ccec9238c46678c38ec7bf6c513b5300", "5c10b6fec79d0412bc07f47dd7c4e8cf"),
               (HemiBMP,       "e13ae1a5a94a1f67decb3c566937e4f9", "a0f3540ba9ab3d98fde3c4f19f33379c"),
               (HexCer,        "f264813c716b565bfc53b4ec5ac5a397", "178e0fa1b427b8afaa97e9a0fb2097b7"),
               (LacCer,        "31c7644b4cd115ee3aaf47fa264782a0", "c78fdfc7108932c2a136d0aa980a0b7f"),
               (LCB_P,         "4d2ba51c5a4d7ee56550168839e892c3", "fb7f5e5a21ff68d2b990ddaa2ac27dd3"),
               (LCB,           "06a29df70c7616fcf0b1c1b7ec068d0f", "824a311eab6fa1fd335c78895af0bc3d"),
               (LPA,           "3ba01733278f450490aae45f36380180", "0e89e1616ca91109cd0b2a5f1fa734cf"),
               (LPC,           "39c1c5ba0c6e703db34338a28d3772ff", "61293d164baaed83f55b70150dfe032d"),
               (LPE,           "824fa10defccf8ffab925f92651a46f9", "9e7e29abf97d3f6d802400e5232592d2"),
               (LPG,           "1c8dfb01edad3fde808d4d6b6ce2fd95", "f83b8a885757b3ca6232965bbe11840c"),
               (LPI,           "d6f79216e5a1600974d01cc44b5a2881", "bae8b6aa95d2595c9ff81567dc411757"),
               (LPS,           "d7cc557847c4774ec36ab6a3b07df606", "2ccab988e601b1756230310f15540c1e"),
               (LysoCL,        "d15b496b7f4c6842b0a9708f2570f677", "c81a4cb958164b36e702d21d35c0a5fb"),
               (LysoCPE,       "34a4e144bccade5ca4cb6d17f640cee7", "2d2091969fa4f518a664adebad8738bb"),
               (LysoCPI,       "7d421f20bf193cc73c3a9da4787819c2", "b98a83744adf806f4f92803676030ae8"),
               (LysoHexCer,    "4afe55ae3b73374e5601506d93b3ced6", "79a9911fed6493dfa810b5370fa9b6e3"),
               (LysoSM,        "bff9bb14293485e04fe6a5309448d80d", "7eb80e66724cc95fd4da8634337a2b20"),
               (MG,            "fdc5aeffded0d8d6b494315d652dee9d", "a72dadb2365c80ea836aab1347b94c3f"),
               (MGDG,          "3ac33c84358a011281645ef85cc291e2", "d8952e957a44532b2d46a4c798ff0535"),
               (MIP2C,         "143edab66a6f3f9d05acdf3b411d72bc", "8f23c164cf2bc31279e8dbcd31e85077"),
               (MIPC,          "0a0045c64d107ee9699de59633eacb8a", "dd23fd98630a595484a6c24827bb79c3"),
               (MMPE,          "17477cc4edf776b0dc4d53c580c68723", "28dffceaaeea8d04ce6ad843b7e0f7f9"),
               (N_Acyl_PE,     "18fcbfdd89c189fd1d44cd7f8625f895", "93654e97567573d9bf4ed53468ce670e"),
               (N_Acyl_PS,     "c992fbb692b98b31a7f90e1880d3aa77", "53aec87f38afd3c0e75f1ad94da3af37"),
               (PA,            "12ff13248cdd99c3cd901c49486ed012", "f242807993777afb38d51e231aece732"),
               (PC,            "324408c9fd8008952c68d4eeef68bf72", "ac6f89b334ad2a604a414bb181e7a4f8"),
               (PE,            "ae98563b06e849aacee9126a58fbe2f0", "f4a5392929ccf38cbbadec8e74905334"),
               (PG,            "129c904faa526ac02e4ad62c66638f71", "d6518c7764afa25992158ae5dafa2d73"),
               (PI,            "a99fbbaf84ef71dd1928e7b8d86afca0", "f51a1db8b54513673a6984a070ef21b0"),
               (PS,            "debaa560107fea3fe8756a996c9d6d78", "5e4d3a20bc47c364fa04eb4271783456"),
               (SM,            "ef05a57b9f52024d3e1c650665517ae6", "e2b393b592254f5501e9333d228c8822"),
               (Sulfatide,     "8b6e01cbaac0df28e6d1d21c9e1ebb58", "c788abd95414fc52429ba926e2554ce5"),
               (Taurine,       "37180d232a0ea95905760509a7d9321a", "06ddff720597fe6e3b482decc2d97258"),
               (TG,            "15704abc9bd04520218e377bada23c90", "349f50d95b5fc0317223ea37f557f74b"),
               ]


@pytest.mark.parametrize("test_input,expected_pos,expected_neg", test_params)
def test_generated_libraries(test_input, expected_pos, expected_neg):
    pos = checksum_for_generated_library(test_input, "pos")
    neg = checksum_for_generated_library(test_input, "neg")
    # print '({:<14} "{}", "{}"), '.format(test_input.__name__ + ",", pos, neg)
    assert pos == expected_pos
    assert neg == expected_neg


