use crate::alphabet::DUMMY_CODE;
use ksw2::{ksw_extz, ksw_extz_t, ksw_gg2, KSW_EZ_RIGHT, KSW_EZ_SCORE_ONLY};
use std::ptr;
use structopt::StructOpt;

#[derive(Clone, StructOpt)]
pub struct AlignmentConfig {
    #[structopt(long = "ma", default_value)]
    pub match_score: i8,
    #[structopt(long = "mp", default_value)]
    pub mismatch_score: i8,
    #[structopt(long = "go", default_value)]
    pub gap_open_penalty: i8,
    #[structopt(long = "ge", default_value)]
    pub gap_extend_penalty: i8,
}

impl Default for AlignmentConfig {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: -4,
            gap_open_penalty: 5,
            gap_extend_penalty: 3,
        }
    }
}

pub struct Aligner {
    config: AlignmentConfig,
    alphabet_size: i8,
    score_matrix: Vec<i8>,
}

impl Aligner {
    pub fn new(config: AlignmentConfig) -> Self {
        let n = DUMMY_CODE as usize + 1;
        let mut score_matrix = vec![0; n * n];
        for i in 1..n - 1 {
            for j in 1..n - 1 {
                score_matrix[i + j * n] = if i == j {
                    config.match_score
                } else {
                    config.mismatch_score
                };
            }
        }

        Self {
            config,
            alphabet_size: n as i8,
            score_matrix,
        }
    }

    pub fn config(&self) -> &AlignmentConfig {
        &self.config
    }

    pub fn global_align(&self, query: &[u8], target: &[u8]) -> i32 {
        unsafe {
            ksw_gg2(
                ptr::null_mut(),
                query.len() as i32,
                query.as_ptr(),
                target.len() as i32,
                target.as_ptr(),
                self.alphabet_size,
                self.score_matrix.as_ptr(),
                self.config.gap_open_penalty,
                self.config.gap_extend_penalty,
                -1,
                ptr::null_mut(),
                ptr::null_mut(),
                ptr::null_mut(),
            )
        }
    }

    pub fn extension_align(&self, query: &[u8], target: &[u8]) -> i32 {
        let mut ez: ksw_extz_t = unsafe { std::mem::zeroed() };
        unsafe {
            ksw_extz(
                ptr::null_mut(),
                query.len() as i32,
                query.as_ptr(),
                target.len() as i32,
                target.as_ptr(),
                self.alphabet_size,
                self.score_matrix.as_ptr(),
                self.config.gap_open_penalty,
                self.config.gap_extend_penalty,
                -1,
                -1,
                (KSW_EZ_SCORE_ONLY | KSW_EZ_RIGHT) as i32,
                &mut ez,
            );
        }
        ez.mqe
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet;

    const MA: i32 = 2;
    const MP: i32 = -4;
    const GO: i32 = 5;
    const GE: i32 = 3;
    const CONFIG: AlignmentConfig = AlignmentConfig {
        match_score: MA as i8,
        mismatch_score: MP as i8,
        gap_open_penalty: GO as i8,
        gap_extend_penalty: GE as i8,
    };

    #[test]
    fn global_align() {
        let aligner = Aligner::new(CONFIG);
        let target = alphabet::encode_seq(b"atcgggatatatggagagcttagag");

        let query1 = alphabet::encode_seq(b"atcgggatatatggagagcttagag");
        assert_eq!(
            aligner.global_align(&query1, &target),
            MA * query1.len() as i32
        );

        let query2 = alphabet::encode_seq(b"atcgggatata");
        assert_eq!(
            aligner.global_align(&query2, &target),
            MA * query2.len() as i32 - GO - GE * (target.len() - query2.len()) as i32
        );
    }

    #[test]
    fn extension_align() {
        let aligner = Aligner::new(CONFIG);
        let target = alphabet::encode_seq(b"atcgggatatatggagagcttagag");
        let query = alphabet::encode_seq(b"atcgggatata");
        assert_eq!(
            aligner.extension_align(&query, &target),
            MA * query.len() as i32,
        );
    }
}
